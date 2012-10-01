package intra;

import java.awt.image.BufferedImage;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Scanner;

import javax.imageio.ImageIO;

import org.ejml.alg.dense.linsol.LinearSolver;
import org.ejml.alg.dense.linsol.svd.SolvePseudoInverseSvd;
import org.ejml.data.DenseMatrix64F;

public class Intra2 {

	public static final int MODES=10;

	public static final int DISK_BLOCK_SIZE=4096;

	public static final int STEPS=30;

	public static final boolean USE_CG=true;

	public static final int[] MODE_COLORS={
		0xFF000000,
		0xFFFFFFFF,
		0xFFFF0000,
		0xFFFFBF00,
		0xFF80FF00,
		0xFF00FF3F,
		0xFF00FFFF,
		0xFF0040FF,
		0xFF7F00FF,
		0xFFFF00C0,
	};

	protected final int B_SZ;

	protected final String DATA_FOLDER;

	protected final int[] INDEX;

	class ModeData {

		protected long numBlocks;

		protected double weightTotal;

		protected double[] mean=new double[2*B_SZ*2*B_SZ];

		protected double[][] covariance=new double[2*B_SZ*2*B_SZ][2*B_SZ*2*B_SZ];

		protected double[][] beta_1=new double[3*B_SZ*B_SZ][B_SZ*B_SZ];

		protected double[] beta_0=new double[B_SZ*B_SZ];

		protected double[] mse=new double[B_SZ*B_SZ];

		protected void addBlock(double _weight,int[] _data) {
			update(_weight,_data);
			numBlocks++;
		}

		protected void removeBlock(double _weight,int[] _data) {
			update(-_weight,_data);
			numBlocks--;
		}

		private void update(double _weight,int[] _data) {
			double[] delta=new double[2*B_SZ*2*B_SZ];
			weightTotal+=_weight;
			// online update of the mean
			for (int i=0;i<2*B_SZ*2*B_SZ;i++) {
				delta[i]=_data[INDEX[i]]-mean[i];
				mean[i]+=delta[i]*_weight/weightTotal;
			}
			// online update of the covariance
			for (int j=0;j<2*B_SZ*2*B_SZ;j++) {
				for (int i=0;i<2*B_SZ*2*B_SZ;i++) {
					covariance[j][i]+=delta[j]*delta[i]*_weight*(weightTotal-_weight)/weightTotal;
				}
			}
		}

		protected void computeBetas() {
			DenseMatrix64F xtx=new DenseMatrix64F(3*B_SZ*B_SZ,3*B_SZ*B_SZ);
			DenseMatrix64F xty=new DenseMatrix64F(3*B_SZ*B_SZ,B_SZ*B_SZ);
			DenseMatrix64F b=new DenseMatrix64F(3*B_SZ*B_SZ,B_SZ*B_SZ);

			// extract X^T*X and X^T*Y
			for (int j=0;j<3*B_SZ*B_SZ;j++) {
				for (int i=0;i<3*B_SZ*B_SZ;i++) {
					xtx.set(j,i,covariance[j][i]/weightTotal);
				}
				for (int i=0;i<B_SZ*B_SZ;i++) {
					xty.set(j,i,covariance[j][3*B_SZ*B_SZ+i]/weightTotal);
				}
			}

			// compute b = (X^T*X)^-1 * X^T*Y 
			LinearSolver<DenseMatrix64F> solver=new SolvePseudoInverseSvd();
			solver.setA(xtx);
			solver.solve(xty,b);

			// extract beta_1
			// compute beta_0 = Y - X * beta_1
			// compute MSE = Y^T*Y - Y^T*X * beta_1
			for (int i=0;i<B_SZ*B_SZ;i++) {
				beta_0[i]=mean[3*B_SZ*B_SZ+i];
				mse[i]=covariance[3*B_SZ*B_SZ+i][3*B_SZ*B_SZ+i]/weightTotal;
				for (int j=0;j<3*B_SZ*B_SZ;j++) {
					beta_1[j][i]=b.get(j,i);
					beta_0[i]-=beta_1[j][i]*mean[j];
					mse[i]-=covariance[3*B_SZ*B_SZ+i][j]/weightTotal*beta_1[j][i];
				}
			}
		}

		protected double msePerCoeff(double[] _mse) {
			double total_mse=0;
			for (int j=0;j<B_SZ*B_SZ;j++) {
				total_mse+=_mse[j];
			}
			return(total_mse/(B_SZ*B_SZ));
		}

		protected double cgPerCoeff(double[] _mse) {
			double total_cg=0;
			for (int j=0;j<B_SZ*B_SZ;j++) {
				total_cg+=-10*Math.log10(_mse[j]/(covariance[3*B_SZ*B_SZ+j][3*B_SZ*B_SZ+j]/weightTotal));
			}
			return(total_cg/(B_SZ*B_SZ));
		}

		protected void printStats() {
			System.out.println("  "+numBlocks+"\t"+msePerCoeff(mse)+"\t"+cgPerCoeff(mse));
		}

		protected double predError(int[] _data,int _y) {
			double pred=beta_0[_y];
			for (int i=0;i<3*B_SZ*B_SZ;i++) {
				pred+=_data[INDEX[i]]*beta_1[i][_y];
			}
			return(Math.abs(_data[INDEX[3*B_SZ*B_SZ+_y]]-pred));
		}

		protected void mseUpdateHelper(int[] _data,int _weight,double[] _mse) {
			for (int j=0;j<B_SZ*B_SZ;j++) {
				double e=predError(_data,j);
				double se=e*e;
				double wt=weightTotal+_weight; 
				double delta=se-mse[j];
				_mse[j]=mse[j]+delta*_weight/wt;
			}
		}

		protected void mseUpdate(int[] _data,int _weight) {
			mseUpdateHelper(_data,_weight,mse);
		}

		protected double deltaBits(int[] _data,int _weight) {
			double old_cg=cgPerCoeff(mse);
			double[] mse=new double[B_SZ*B_SZ];
			mseUpdateHelper(_data,_weight,mse);
			update(_weight,_data);
			double cg=cgPerCoeff(mse);
			update(-_weight,_data);
			return(cg-old_cg);
		}

	};

	public Intra2(int _blockSize,String _dataFolder) {
		B_SZ=_blockSize;
		DATA_FOLDER=_dataFolder;
		INDEX=new int[4*B_SZ*B_SZ];
		for (int i=0,j=0;j<B_SZ;j++) {
			for (int k=0;k<B_SZ;k++) {
				INDEX[i+0*B_SZ*B_SZ]=2*B_SZ*j+k;
				INDEX[i+1*B_SZ*B_SZ]=2*B_SZ*j+(B_SZ+k);
				INDEX[i+2*B_SZ*B_SZ]=2*B_SZ*(B_SZ+j)+k;
				INDEX[i+3*B_SZ*B_SZ]=2*B_SZ*(B_SZ+j)+(B_SZ+k);
				i++;
			}
		}
	}

	protected File[] getFiles() throws IOException {
		File folder=new File(DATA_FOLDER);
		File[] coeffFiles=folder.listFiles(new FilenameFilter() {
			public boolean accept(File _file,String _filename) {
				return(_filename.endsWith(".coeffs"));
			}
		});
		if (coeffFiles==null) {
			return(null);
		}
		Arrays.sort(coeffFiles);
		File[] binFiles=new File[coeffFiles.length];
		for (int i=0;i<coeffFiles.length;i++) {
			binFiles[i]=new File(coeffFiles[i].getPath()+".bin");
			if (!binFiles[i].exists()) {
				System.out.println("Converting "+coeffFiles[i].getPath()+" to "+binFiles[i].getPath());
				Scanner s=new Scanner(new BufferedInputStream(new FileInputStream(coeffFiles[i]),DISK_BLOCK_SIZE));
				DataOutputStream dos=new DataOutputStream(new BufferedOutputStream(new FileOutputStream(binFiles[i]),DISK_BLOCK_SIZE));
				while (s.hasNextShort()) {
					dos.writeShort(s.nextShort());
				}
				s.close();
				dos.close();
			}
		}
		return(binFiles);
	}

	protected ModeData[] loadData(File[] _files) throws IOException {
		// initialize the modes
		ModeData modeData[]=new ModeData[MODES];
		for (int i=0;i<MODES;i++) {
			modeData[i]=new ModeData();
		}

		for (File file : _files) {
			System.out.println("Loading "+file.getPath());
			DataInputStream dis=new DataInputStream(new BufferedInputStream(new FileInputStream(file),DISK_BLOCK_SIZE));
			DataOutputStream dos=new DataOutputStream(new BufferedOutputStream(new FileOutputStream(file.getPath()+".step00.mode"),DISK_BLOCK_SIZE));
			// 3 color planes, YUV
			for (int pli=0;pli<1;pli++) {
				// read the plane size
				int nx=dis.readShort();
				int ny=dis.readShort();

				int[] rgb=new int[nx*ny];
				for (int block=0;block<nx*ny;block++) {
					int mode=dis.readShort();
					int weight=dis.readShort();
					//System.out.println("mode="+mode+" weight="+weight);
					dos.writeShort(mode);
					//dos.writeShort(weight);
					int[] data=new int[2*B_SZ*2*B_SZ];
					for (int i=0;i<2*B_SZ*2*B_SZ;i++) {
						data[i]=dis.readShort();
					}
					if (weight>0) {
						modeData[mode].addBlock(mode==0?1:weight,data);
						//modeData[mode].addBlock(weight,data);
					}
					rgb[block]=MODE_COLORS[mode];
				}

				BufferedImage bi=new BufferedImage(nx,ny,BufferedImage.TYPE_INT_ARGB);
				bi.setRGB(0,0,nx,ny,rgb,0,nx);
				ImageIO.write(bi,"PNG",new File(file.getPath()+".step00.png"));
			}
			dos.close();
			dis.close();
		}

		return(modeData);
	}

	protected void processData(int step,ModeData[] _modeData,File[] _files) throws IOException {
		for (File file : _files) {
			System.out.println("Processing "+file.getPath());
			DataInputStream dis=new DataInputStream(new BufferedInputStream(new FileInputStream(file),DISK_BLOCK_SIZE));
			DataInputStream dis2=new DataInputStream(new BufferedInputStream(new FileInputStream(file.getPath()+".step"+(step-1<10?"0"+(step-1):step-1)+".mode"),DISK_BLOCK_SIZE));
			DataOutputStream dos=new DataOutputStream(new BufferedOutputStream(new FileOutputStream(file.getPath()+".step"+(step<10?"0"+step:step)+".mode"),DISK_BLOCK_SIZE));

			// read the plane size
			int nx=dis.readShort();
			int ny=dis.readShort();

			int[] rgb=new int[nx*ny];
			for (int block=0;block<nx*ny;block++) {
				// load the data
				int mode=dis.readShort();
				int weight=dis.readShort();
				int[] data=new int[2*B_SZ*2*B_SZ];
				for (int i=0;i<2*B_SZ*2*B_SZ;i++) {
					data[i]=dis.readShort();
				}
				int lastMode=dis2.readShort();
				if (weight>0) {
					// compute error
					double[] error=new double[MODES];
					double[] cg=new double[MODES];
					for (int i=0;i<MODES;i++) {
						for (int j=0;j<B_SZ*B_SZ;j++) {
							error[i]+=_modeData[i].predError(data,j);
						}
						if (USE_CG) {
							cg[i]=_modeData[i].deltaBits(data,(lastMode==i?-1:1)*(i==0?1:weight));
						}
					}
					double best=-Double.MAX_VALUE;
					double nextBest=-Double.MAX_VALUE;
					for (int j=0;j<MODES;j++) {
						// lower SATD is better
						double metric=error[lastMode]-error[j];
						if (USE_CG) {
							// greater CG is better
							metric=lastMode==j?0:cg[lastMode]+cg[j];
						}
						if (metric>best) {
							nextBest=best;
							best=metric;
							mode=j;
						}
						else {
							if (metric>nextBest) {
								nextBest=metric;
							}
						}
					}
					if (mode!=lastMode) {
						_modeData[lastMode].mseUpdate(data,lastMode==0?-1:-weight);
						_modeData[mode].mseUpdate(data,mode==0?1:weight);
						_modeData[lastMode].removeBlock(lastMode==0?1:weight,data);
						_modeData[mode].addBlock(mode==0?1:weight,data);
					}
				}
				dos.writeShort(mode);
				rgb[block]=MODE_COLORS[mode];
			}

			BufferedImage bi=new BufferedImage(nx,ny,BufferedImage.TYPE_INT_ARGB);
			bi.setRGB(0,0,nx,ny,rgb,0,nx);
			ImageIO.write(bi,"PNG",new File(file.getPath()+".step"+(step<10?"0"+step:step)+".png"));

			dis.close();
			dis2.close();
			dos.close();
		}
	}

	public void run() throws Exception {
		File[] files=getFiles();
		if (files==null) {
			System.out.println("no data files in "+DATA_FOLDER);
			return;
		}

		// load the blocks
		ModeData[] modeData=loadData(files);

		long start=System.currentTimeMillis();
		for (int k=1;k<=STEPS;k++) {
			// compute betas and MSE
			for (int i=0;i<MODES;i++) {
				modeData[i].computeBetas();
				modeData[i].printStats();
			}

			// reclassify blocks
			processData(k,modeData,files);

			long now=System.currentTimeMillis();
			System.out.println("Step "+k+" took "+((now-start)/1000.0)+"s");
		}
	}

	public static void main(String[] _args) throws Exception {
		Intra2 intra=new Intra2(4,"data");
		intra.run();
	}

};
