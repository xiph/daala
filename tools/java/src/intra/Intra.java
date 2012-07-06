package intra;

import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;

import javax.imageio.ImageIO;

import org.ejml.alg.dense.linsol.LinearSolver;
import org.ejml.alg.dense.linsol.svd.SolvePseudoInverseSvd;
import org.ejml.data.DenseMatrix64F;

public class Intra {

	public static final int MODES=10;

	public static final int B_SZ=4;

	public static final int MAX_MULTS=4*B_SZ*B_SZ;

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

	public static final int[] INDEX;

	static {
		// index the coefficients in block order
		INDEX=new int[4*B_SZ*B_SZ];
		for (int tmp=0,j=0;j<B_SZ;j++) {
			for (int i=0;i<B_SZ;i++) {
				INDEX[tmp+0*B_SZ*B_SZ]=2*B_SZ*j+i;
				INDEX[tmp+1*B_SZ*B_SZ]=2*B_SZ*j+(B_SZ+i);
				INDEX[tmp+2*B_SZ*B_SZ]=2*B_SZ*(B_SZ+j)+i;
				INDEX[tmp+3*B_SZ*B_SZ]=2*B_SZ*(B_SZ+j)+(B_SZ+i);
				tmp++;
			}
		}


	}
	static class ModeData {

		protected long numBlocks;

		protected double weightTotal;

		protected double[] mean=new double[2*B_SZ*2*B_SZ];

		protected double[][] covariance=new double[2*B_SZ*2*B_SZ][2*B_SZ*2*B_SZ];

		protected double[] scale=new double[2*B_SZ*2*B_SZ];

		protected int[] numMults=new int[B_SZ*B_SZ];

		protected boolean[][] mult=new boolean[3*B_SZ*B_SZ][B_SZ*B_SZ];

		protected double[][] mseMasked=new double[3*B_SZ*B_SZ][B_SZ*B_SZ];

		protected double[] mse=new double[B_SZ*B_SZ];

		protected double mseTotal;

		protected double[][] beta=new double[3*B_SZ*B_SZ][B_SZ*B_SZ];

	};

	protected static ModeData[] modeData=new ModeData[MODES];

	protected static double mse(int _mode,int _y) {
		return(mse(_mode,_y,false));
	}

	protected static LinearSolver<DenseMatrix64F> solver=new SolvePseudoInverseSvd();

	protected static DenseMatrix64F[] xtx=new DenseMatrix64F[3*B_SZ*B_SZ];

	protected static DenseMatrix64F[] xty=new DenseMatrix64F[3*B_SZ*B_SZ];

	protected static DenseMatrix64F[] beta=new DenseMatrix64F[3*B_SZ*B_SZ];

	protected static double mse(int _mode,int _y,boolean _saveBeta) {
		ModeData md=modeData[_mode];
		double yty=md.covariance[3*B_SZ*B_SZ+_y][3*B_SZ*B_SZ+_y];
		double ytxb=0;
		int xsize=md.numMults[_y];
		if (xsize>0) {
			if (xtx[xsize-1]==null) {
				xty[xsize-1]=new DenseMatrix64F(xsize,1);
				beta[xsize-1]=new DenseMatrix64F(xsize,1);
			}
			xtx[xsize-1]=new DenseMatrix64F(xsize,xsize);
			for (int ji=0,j=0;j<3*B_SZ*B_SZ;j++) {
				if (md.mult[j][_y]) {
					for (int ii=0,i=0;i<3*B_SZ*B_SZ;i++) {
						if (md.mult[i][_y]) {
							xtx[xsize-1].set(ji,ii,md.covariance[j][i]);
							ii++;
						}
					}
					xty[xsize-1].set(ji,0,md.covariance[j][3*B_SZ*B_SZ+_y]);
					ji++;
				}
			}

			solver.setA(xtx[xsize-1]);
			solver.solve(xty[xsize-1],beta[xsize-1]);
			// compute y^T * x * beta
			for (int ii=0,i=0;i<3*B_SZ*B_SZ;i++) {
				if (md.mult[i][_y]) {
					ytxb+=md.covariance[3*B_SZ*B_SZ+_y][i]*beta[xsize-1].get(ii,0);
					ii++;
				}
			}
			// save beta for use with classification
			if (_saveBeta) {
				for (int ii=0,i=0;i<3*B_SZ*B_SZ;i++) {
					if (md.mult[i][_y]) {
						// beta' = Sx^-1 * beta * Sy -> beta = Sx * beta' * Sy^-1
						md.beta[i][_y]=beta[xsize-1].get(ii,0)*md.scale[3*B_SZ*B_SZ+_y]/md.scale[i];
						ii++;
					}
					else {
						md.beta[i][_y]=0;
					}
				}
			}
		}
		return(yty-ytxb);
	}

	// compute the MSE for a given y coefficient
	protected static void mseUpdate(int _mode,int _y) {
		ModeData md=modeData[_mode];
		for (int i=0;i<3*B_SZ*B_SZ;i++) {
			if (md.mult[i][_y]) {
				md.mult[i][_y]=false;
				md.numMults[_y]--;
				md.mseMasked[i][_y]=mse(_mode,_y);
				md.mult[i][_y]=true;
				md.numMults[_y]++;
			 }
		 }
	}

	protected static int mseDrop(int _mode) {
		ModeData md=modeData[_mode];
		double delta=Double.MAX_VALUE;
		int y=-1;
		int c=-1;
		for (int j=0;j<3*B_SZ*B_SZ;j++) {
			for (int i=0;i<B_SZ*B_SZ;i++) {
				if (md.mult[j][i]) {
					if (md.mseMasked[j][i]-md.mse[i]<delta) {
						delta=md.mseMasked[j][i]-md.mse[i];
						y=i;
						c=j;
					}
				}
			}
		}
		md.mult[c][y]=false;
		return(y);
	}

	protected static void addBlock(int _mode,double _weight,String[] _data) {
		double[] delta=new double[2*B_SZ*2*B_SZ];
		modeData[_mode].numBlocks++;
		modeData[_mode].weightTotal+=_weight;
		// online update of the mean
		for (int i=0;i<2*B_SZ*2*B_SZ;i++) {
			delta[i]=Integer.parseInt(_data[INDEX[i]+1])-modeData[_mode].mean[i];
			modeData[_mode].mean[i]+=delta[i]*_weight/modeData[_mode].weightTotal;
		}
		// online update of the covariance
		for (int j=0;j<2*B_SZ*2*B_SZ;j++) {
			for (int i=0;i<2*B_SZ*2*B_SZ;i++) {
				modeData[_mode].covariance[j][i]+=delta[j]*delta[i]*_weight*(modeData[_mode].weightTotal-_weight)/modeData[_mode].weightTotal;
			}
		}
	}

	public static void main(String[] _args) throws IOException {
		// initialize the modes
		for (int i=0;i<MODES;i++) {
			modeData[i]=new ModeData();
		}

		// for each pass of k-means load the data
		File folder=new File("data");
		File[] inputFiles=folder.listFiles(new FilenameFilter() {
			public boolean accept(File _file,String _filename) {
				return(_filename.endsWith(".coeffs"));
			}
		});
		for (File file : inputFiles) {
			System.out.println(file.getPath());
			BufferedReader br=new BufferedReader(new FileReader(file));
			// 3 color planes, YUV
			for (int pli=0;pli<1;pli++) {
				// read the plane size
				String[] data=br.readLine().split(" ");
				int nx=Integer.parseInt(data[0]);
				int ny=Integer.parseInt(data[1]);

				int[] rgb=new int[nx*ny];
				for (int block=0;block<nx*ny;block++) {
					data=br.readLine().split(" ");
					int mode=Integer.parseInt(data[0]);
					addBlock(mode,1,data);
					rgb[block]=MODE_COLORS[mode];
				}

				BufferedImage bi=new BufferedImage(nx,ny,BufferedImage.TYPE_INT_ARGB);
				bi.setRGB(0,0,nx,ny,rgb,0,nx);
				ImageIO.write(bi,"PNG",new File(file.getPath()+".step0.png"));
			}
		}

		long last=System.currentTimeMillis();
		// 16*44=704=768-64
		for (int step=1;step<=44;step++) {
			long modeStart=System.currentTimeMillis();
			// for each mode
			for (int mode=0;mode<MODES;mode++) {
				ModeData md=modeData[mode];
				for (int i=0;i<2*B_SZ*2*B_SZ;i++) {
					md.scale[i]=Math.sqrt(md.covariance[i][i]);
				}
				// compute C' = Sx * C * Sx
				for (int j=0;j<2*B_SZ*2*B_SZ;j++) {
					for (int i=0;i<2*B_SZ*2*B_SZ;i++) {
						md.covariance[j][i]/=md.scale[j]*md.scale[i];
					}
				}
				// compute MSE based on completely unmasked
				md.mseTotal=0;
				for (int y=0;y<B_SZ*B_SZ;y++) {
					md.numMults[y]=3*B_SZ*B_SZ;
					for (int i=0;i<3*B_SZ*B_SZ;i++) {
						md.mult[i][y]=true;
					}
					md.mse[y]=mse(mode,y,true);
					md.mseTotal+=md.mse[y];
					mseUpdate(mode,y);
				}

				System.out.println("before="+(System.currentTimeMillis()-modeStart));
				double rmseStart=Math.sqrt(md.mseTotal/(B_SZ*B_SZ));
				for (int drops=0;drops<step*16;drops++) {
					int y=mseDrop(mode);
					md.mseTotal-=md.mse[y];
					md.mse[y]=mse(mode,y,true);
					md.mseTotal+=md.mse[y];
					md.numMults[y]--;
					mseUpdate(mode,y);
				}
				System.out.println("after="+(System.currentTimeMillis()-modeStart));
				double rmseEnd=Math.sqrt(md.mseTotal/(B_SZ*B_SZ));
				long modeEnd=System.currentTimeMillis();
				System.out.println(mode+" "+rmseEnd+" "+(rmseEnd-rmseStart)+" took "+(modeEnd-modeStart)/1000.0);
				modeStart=modeEnd;
				for (int j=0;j<B_SZ;j++) {
					String prefix="    ";
					for (int i=0;i<B_SZ;i++) {
						System.out.print(prefix+(md.numMults[j*B_SZ+i]));
						prefix=" ";
					}
					prefix="    ";
					for (int i=0;i<B_SZ;i++) {
						System.out.print(prefix+md.mse[j*B_SZ+i]);
						prefix=" ";
					}
					System.out.println();
				}
			}

			// reset the mode data
			for (ModeData md : modeData) {
				md.numBlocks=0;
				md.weightTotal=0;
				for (int j=0;j<2*B_SZ*2*B_SZ;j++) {
					md.mean[j]=0;
					for (int i=0;i<2*B_SZ*2*B_SZ;i++) {
						md.covariance[j][i]=0;
					}
				}
			}

			// re-classify the blocks based on the computed beta's
			double satdTotal=0;
			double blocksTotal=0;
			for (File file : inputFiles) {
				System.out.println(file.getPath());
				BufferedReader br=new BufferedReader(new FileReader(file));
				// 3 color planes, YUV
				for (int pli=0;pli<1;pli++) {
					// read the plane size
					String[] data=br.readLine().split(" ");
					int nx=Integer.parseInt(data[0]);
					int ny=Integer.parseInt(data[1]);

					int[] rgb=new int[nx*ny];
					for (int block=0;block<nx*ny;block++) {
						data=br.readLine().split(" ");
						double best=Double.MAX_VALUE;
						int mode=-1;
						double nextBest=0;
						for (int j=0;j<MODES;j++) {
							double satd=0;
							// compute |y-x*B| the SATD (or L1-norm)
							for (int y=0;y<B_SZ*B_SZ;y++) {
								double pred=0;
								for (int i=0;i<3*B_SZ*B_SZ;i++) {
									pred+=Integer.parseInt(data[INDEX[i]+1])*modeData[j].beta[i][y];
								}
								satd+=Math.abs(Integer.parseInt(data[INDEX[3*B_SZ*B_SZ+y]+1])-pred);
							}
							if (satd<best) {
								nextBest=best;
								best=satd;
								mode=j;
							}
						}
						satdTotal+=best;
						addBlock(mode,mode==0?1:nextBest-best,data);
						rgb[block]=MODE_COLORS[mode];
					}
					blocksTotal+=nx*ny;

					BufferedImage bi=new BufferedImage(nx,ny,BufferedImage.TYPE_INT_ARGB);
					bi.setRGB(0,0,nx,ny,rgb,0,nx);
					ImageIO.write(bi,"PNG",new File(file.getPath()+".step"+(step<10?"0"+step:step)+".png"));
				}
			}
			long now=System.currentTimeMillis();
			System.out.println("Step "+step+" took "+(now-last)/1000.0);
			System.out.println("Average SATD: "+(satdTotal/blocksTotal));
		}
	}

};
