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

	public static final int LEN=(2*B_SZ)*(2*B_SZ);

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

	static class ModeData {

		protected long numBlocks;

		protected double[] mean=new double[LEN];

		protected double[][] covariance=new double[LEN][LEN];

		protected int[] numMasked=new int[B_SZ*B_SZ];

		protected boolean[][] mask=new boolean[3*B_SZ*B_SZ][B_SZ*B_SZ];

		protected double[][] rmseMasked=new double[3*B_SZ*B_SZ][B_SZ*B_SZ];

		protected double[] rmse=new double[B_SZ*B_SZ];

		protected double rmseTotal;

	};

	protected static ModeData[] modeData=new ModeData[MODES];

	protected static double rmse(int _mode,int _y) {
		ModeData md=modeData[_mode];
		double yty=md.covariance[3*B_SZ*B_SZ+_y][3*B_SZ*B_SZ+_y];
		double ytxb=0;
		int xsize=3*B_SZ*B_SZ-md.numMasked[_y];
		if (xsize>0) {
			DenseMatrix64F xtx=new DenseMatrix64F(xsize,xsize);
			DenseMatrix64F xty=new DenseMatrix64F(xsize,1);
			for (int ji=0,j=0;j<3*B_SZ*B_SZ;j++) {
				if (!md.mask[j][_y]) {
					for (int ii=0,i=0;i<3*B_SZ*B_SZ;i++) {
						if (!md.mask[i][_y]) {
							xtx.set(ji,ii,md.covariance[j][i]);
							ii++;
						}
					}
					xty.set(ji,0,md.covariance[j][3*B_SZ*B_SZ+_y]);
					ji++;
				}
			}
			LinearSolver<DenseMatrix64F> solver=new SolvePseudoInverseSvd();
			DenseMatrix64F beta=new DenseMatrix64F(xsize,1);
			solver.setA(xtx);
			solver.solve(xty,beta);
			for (int ii=0,i=0;i<3*B_SZ*B_SZ;i++) {
				if (!md.mask[i][_y]) {
					ytxb+=md.covariance[3*B_SZ*B_SZ+_y][i]*beta.get(ii,0);
					ii++;
				}
			}
		}
		return((yty-ytxb)/md.numBlocks);
	}

	// compute the RMSE for a given y coefficient
	protected static void rmseUpdate(int _mode,int _y) {
		ModeData md=modeData[_mode];
		for (int i=0;i<3*B_SZ*B_SZ;i++) {
			if (!md.mask[i][_y]) {
				md.mask[i][_y]=true;
				md.numMasked[_y]++;
				md.rmseMasked[i][_y]=rmse(_mode,_y);
				md.mask[i][_y]=false;
				md.numMasked[_y]--;
			 }
		 }
	}

	protected static void rmsePrint(int _mode,int _mults) {
		ModeData md=modeData[_mode];
		System.out.print("mults["+_mode+"]="+_mults+" rmse["+_mode+"]="+md.rmseTotal);
		for (int y=0;y<B_SZ*B_SZ;y++) {
			System.out.print(" "+md.rmse[y]);
		}
		System.out.println();
	}

	protected static int rmseDrop(int _mode) {
		ModeData md=modeData[_mode];
		double delta=Double.MAX_VALUE;
		int y=-1;
		int c=-1;
		for (int j=0;j<3*B_SZ*B_SZ;j++) {
			for (int i=0;i<B_SZ*B_SZ;i++) {
				if (!md.mask[j][i]) {
					if (md.rmseMasked[j][i]-md.rmse[i]<0) {
						System.out.println("y="+i+" c="+j+" old rmse="+md.rmse[i]+" new rmse="+md.rmseMasked[j][i]);
					}
					if (md.rmseMasked[j][i]-md.rmse[i]<delta) {
						delta=md.rmseMasked[j][i]-md.rmse[i];
						y=i;
						c=j;
					}
				}
			}
		}
		//System.out.println("delta="+delta);
		md.mask[c][y]=true;
		return(y);
	}

	public static void main(String[] _args) throws IOException {
		for (int i=0;i<MODES;i++) {
			modeData[i]=new ModeData();
		}

		// index the coefficients in block order
		int[] index=new int[4*B_SZ*B_SZ];
		int tmp=0;
		for (int j=0;j<B_SZ;j++) {
			for (int i=0;i<B_SZ;i++) {
				index[tmp+0*B_SZ*B_SZ]=2*B_SZ*j+i;
				index[tmp+1*B_SZ*B_SZ]=2*B_SZ*j+(B_SZ+i);
				index[tmp+2*B_SZ*B_SZ]=2*B_SZ*(B_SZ+j)+i;
				index[tmp+3*B_SZ*B_SZ]=2*B_SZ*(B_SZ+j)+(B_SZ+i);
				tmp++;
			}
		}

		// for each pass of k-means load the data
		File folder=new File("data");
		File[] inputFiles=folder.listFiles(new FilenameFilter() {
			public boolean accept(File _file,String _filename) {
				return(_filename.endsWith(".coeffs"));
			}
		});
		for (File file : inputFiles) {
			System.out.println("Processing "+file.getPath());
			BufferedReader br=new BufferedReader(new FileReader(file));
			// 3 color planes, YUV
			for (int pli=0;pli<1;pli++) {
				// read the plane size
				String[] data=br.readLine().split(" ");
				int nx=Integer.parseInt(data[0]);
				int ny=Integer.parseInt(data[1]);

				int[] rgb=new int[nx*ny];
				double[] delta=new double[LEN];
				for (int block=0;block<nx*ny;block++) {
					data=br.readLine().split(" ");
					int mode=Integer.parseInt(data[0]);
					rgb[block]=MODE_COLORS[mode];
					modeData[mode].numBlocks++;
					// online update of the mean
					for (int i=0;i<LEN;i++) {
						delta[i]=Integer.parseInt(data[index[i]+1])-modeData[mode].mean[i];
						modeData[mode].mean[i]+=delta[i]/modeData[mode].numBlocks;
					}
					// online update of the covariance
					for (int j=0;j<LEN;j++) {
						for (int i=0;i<LEN;i++) {
							modeData[mode].covariance[j][i]+=delta[j]*delta[i]*(modeData[mode].numBlocks-1)/modeData[mode].numBlocks;
						}
					}
				}

				BufferedImage bi=new BufferedImage(nx,ny,BufferedImage.TYPE_INT_ARGB);
				bi.setRGB(0,0,nx,ny,rgb,0,nx);
				ImageIO.write(bi,"PNG",new File(file.getPath()+".png"));
			}
		}

		for (int step=1;step<=1;step++) {
			// for each mode
			for (int mode=0;mode<MODES;mode++) {
				ModeData md=modeData[mode];
				// initialize the rmse based on completely unmasked
				md.rmseTotal=0;
				for (int y=0;y<B_SZ*B_SZ;y++) {
					md.rmse[y]=rmse(mode,y);
					md.rmseTotal+=md.rmse[y];
				}
				for (int y=0;y<B_SZ*B_SZ;y++) {
					rmseUpdate(mode,y);
				}
				// print the current rmse
				int mults=(3*B_SZ*B_SZ)*(B_SZ*B_SZ);
				//rmsePrint(mode,mults);
				double rmseStart=Math.sqrt(md.rmseTotal/(B_SZ*B_SZ));
				for (;mults-->MAX_MULTS;) {
					int y=rmseDrop(mode);
					md.rmseTotal-=md.rmse[y];
					md.rmse[y]=rmse(mode,y);
					md.rmseTotal+=md.rmse[y];
					md.numMasked[y]++;
					rmseUpdate(mode,y);
					//rmsePrint(mode,mults);
				}
				double rmseEnd=Math.sqrt(md.rmseTotal/(B_SZ*B_SZ));
				System.out.print(mode+" "+rmseEnd+" "+(rmseEnd-rmseStart));
				for (int i=0;i<B_SZ*B_SZ;i++) {
					System.out.print(" "+(3*B_SZ*B_SZ-md.numMasked[i]));
				}
				System.out.println();
			}
		}
	}

};
