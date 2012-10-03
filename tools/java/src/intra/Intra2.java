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

	public static final int BITS_PER_COEFF=6;

	public static final double MAX_CG=BITS_PER_COEFF*-10*Math.log10(0.5);

	public static final boolean USE_CG=true;

	public static final boolean UPDATE_WEIGHT=false;

	public static double DC_WEIGHT=66.07307086571254;

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
				double cg=-10*Math.log10(_mse[j]/(covariance[3*B_SZ*B_SZ+j][3*B_SZ*B_SZ+j]/weightTotal));
				/*if (cg>=MAX_CG) {
					cg=MAX_CG;
				}*/
				total_cg+=cg;
			}
			return(total_cg/(B_SZ*B_SZ));
		}

		protected double predError(int[] _data,int _y) {
			double pred=beta_0[_y];
			for (int i=0;i<3*B_SZ*B_SZ;i++) {
				pred+=_data[INDEX[i]]*beta_1[i][_y];
			}
			return(Math.abs(_data[INDEX[3*B_SZ*B_SZ+_y]]-pred));
		}

		protected void mseUpdateHelper(int[] _data,double _weight,double[] _mse) {
			for (int j=0;j<B_SZ*B_SZ;j++) {
				double e=predError(_data,j);
				double se=e*e;
				double wt=weightTotal+_weight; 
				double delta=se-mse[j];
				_mse[j]=mse[j]+delta*_weight/wt;
			}
		}

		protected void mseUpdate(int[] _data,double _weight) {
			mseUpdateHelper(_data,_weight,mse);
		}

		protected double deltaBits(int[] _data,double _weight) {
			double old_cg=cgPerCoeff(mse);
			double[] tmp_mse=new double[B_SZ*B_SZ];
			mseUpdateHelper(_data,_weight,tmp_mse);
			update(_weight,_data);
			double cg=cgPerCoeff(tmp_mse);
			update(-_weight,_data);
			return((weightTotal+_weight)*cg-weightTotal*old_cg);
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

		double weight_sum=0;
		long block_sum=0;

		for (File file : _files) {
			System.out.println("Loading "+file.getPath());
			DataInputStream dis=new DataInputStream(new BufferedInputStream(new FileInputStream(file),DISK_BLOCK_SIZE));
			DataOutputStream dos=new DataOutputStream(new BufferedOutputStream(new FileOutputStream(file.getPath()+".step00.mode"),DISK_BLOCK_SIZE));

			// read the plane size
			int nx=dis.readShort();
			int ny=dis.readShort();

			int[] rgb=new int[nx*ny];
			for (int block=0;block<nx*ny;block++) {
				int mode=dis.readShort();
				double weight=dis.readShort();
				//System.out.println("mode="+mode+" weight="+weight);
				int[] data=new int[2*B_SZ*2*B_SZ];
				for (int i=0;i<2*B_SZ*2*B_SZ;i++) {
					data[i]=dis.readShort();
				}
				if (weight>0) {
					weight_sum+=weight;
					block_sum++;
					if (mode==0) {
						weight=DC_WEIGHT;
					}
					modeData[mode].addBlock(weight,data);
				}
				dos.writeShort(mode);
				dos.writeDouble(weight);
				rgb[block]=MODE_COLORS[mode];
			}

			BufferedImage bi=new BufferedImage(nx,ny,BufferedImage.TYPE_INT_ARGB);
			bi.setRGB(0,0,nx,ny,rgb,0,nx);
			ImageIO.write(bi,"PNG",new File(file.getPath()+".step00.png"));

			dos.close();
			dis.close();
		}

		DC_WEIGHT=weight_sum/block_sum;
		System.out.println("DC_WEIGHT="+DC_WEIGHT);

		// correct the DC to use DC_WEIGHT

		return(modeData);
	}

	protected void fitData(ModeData[] _modeData) {
		// compute betas and MSE
		for (int i=0;i<MODES;i++) {
			System.out.println("mode "+i);
			double[] old_mse=new double[B_SZ*B_SZ];
			for (int j=0;j<B_SZ*B_SZ;j++) {
				old_mse[j]=_modeData[i].mse[j];
			}
			_modeData[i].computeBetas();
			for (int j=0;j<B_SZ*B_SZ;j++) {
				System.out.println("  "+j+": "+old_mse[j]+"\t"+_modeData[i].mse[j]+"\t"+(_modeData[i].mse[j]-old_mse[j]));
			}
		}
	}

	protected static final int SPACE=4;

	protected void printStats(ModeData[] _modeData) {
		double mse_sum=0;
		double cg_sum=0;
		double weight=0;
		for (int i=0;i<MODES;i++) {
			double mse=_modeData[i].msePerCoeff(_modeData[i].mse);
			double cg=_modeData[i].cgPerCoeff(_modeData[i].mse);
			System.out.println("  "+i+": "+_modeData[i].numBlocks+"\t"+_modeData[i].weightTotal+"\t"+mse+"\t"+cg);
			mse_sum+=_modeData[i].weightTotal*mse;
			cg_sum+=_modeData[i].weightTotal*cg;
			weight+=_modeData[i].weightTotal;
		}
		System.out.println("Total Weight "+weight);
		System.out.println("Average MSE "+mse_sum/weight);
		System.out.println("Average CG  "+cg_sum/weight);

		// make a "heat map" using the normalized beta_1
		/*int w=SPACE+(2*B_SZ+SPACE)*B_SZ*B_SZ;
		int h=SPACE+(2*B_SZ+SPACE)*MODES;
		int[] rgb=new int[w*h];

		for (int i=0;i<MODES;i++) {
			int yoff=SPACE+i*(2*B_SZ+SPACE);
			int xoff=SPACE;
			for (int j=0;j<B_SZ*B_SZ;j++) {

				rgb[w*(yoff+i)+xoff];
			}
		}*/
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
				if ((block&0x3fff)==0) {
					System.out.println(block);
				}

				// load the data
				int mode=dis.readShort();
				double weight=dis.readShort();
				int[] data=new int[2*B_SZ*2*B_SZ];
				for (int i=0;i<2*B_SZ*2*B_SZ;i++) {
					data[i]=dis.readShort();
				}
				int lastMode=dis2.readShort();
				double lastWeight=dis2.readDouble();
				if (weight>0) {
					if (mode==0) {
						weight=DC_WEIGHT;
						lastWeight=DC_WEIGHT;
					}
					// compute error
					double[] error=new double[MODES];
					double[] cg=new double[MODES];
					for (int i=0;i<MODES;i++) {
						for (int j=0;j<B_SZ*B_SZ;j++) {
							error[i]+=_modeData[i].predError(data,j);
						}
						if (USE_CG) {
							cg[i]=_modeData[i].deltaBits(data,lastMode==i?-weight:weight);
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
					if (UPDATE_WEIGHT) {
						weight=best-nextBest;
					}
					if (USE_CG) {
						_modeData[lastMode].mseUpdate(data,-lastWeight);
						_modeData[mode].mseUpdate(data,weight);
					}
					_modeData[lastMode].removeBlock(lastWeight,data);
					_modeData[mode].addBlock(weight,data);
					/*if (USE_CG) {
						_modeData[lastMode].computeBetas();
						_modeData[mode].computeBetas();
					}*/
				}
				dos.writeShort(mode);
				dos.writeDouble(weight);
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
			System.out.println("No .coeffs files found in "+DATA_FOLDER+" folder.  Enable the PRINT_BLOCKS ifdef and run:");
			System.out.println("  for f in subset1-y4m/*.y4m; do ./init_intra_maps $f && ./init_intra_xform $f 2> $f.coeffs; done");
			System.out.println("Move the *.coeffs files to "+DATA_FOLDER);
			return;
		}

		// load the blocks
		ModeData[] modeData=loadData(files);

		long start=System.currentTimeMillis();
		for (int k=1;k<=STEPS;k++) {

			// update model
			fitData(modeData);

			printStats(modeData);

			// reclassify blocks
			processData(k,modeData,files);

			long now=System.currentTimeMillis();
			System.out.println("Step "+k+" took "+((now-start)/1000.0)+"s");
		}
	}

	public static void main(String[] _args) throws Exception {
		Intra2 intra=new Intra2(4,"data2");
		intra.run();
	}

};
