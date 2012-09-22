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

public class Intra {

	public static final int MODES=10;

	public static final int B_SZ=4;

	public static final boolean DROP_MULTS=false;

	public static final int MAX_MULTS=4*B_SZ*B_SZ;

	public static final boolean USE_CG=false;

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

	public static final int DISK_BLOCK_SIZE=4096;

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

		protected double[][] beta=new double[3*B_SZ*B_SZ][B_SZ*B_SZ];

		protected double satd;

		protected ModeData() {
			for (int y=0;y<B_SZ*B_SZ;y++) {
				numMults[y]=3*B_SZ*B_SZ;
				for (int i=0;i<3*B_SZ*B_SZ;i++) {
					mult[i][y]=true;
				}
			}
		}

		protected void addBlock(double _weight,int[] _data) {
			double[] delta=new double[2*B_SZ*2*B_SZ];
			numBlocks++;
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

		protected void normalize() {
			for (int i=0;i<2*B_SZ*2*B_SZ;i++) {
				scale[i]=Math.sqrt(covariance[i][i]);
			}
			// compute C' = Sx * C * Sx
			for (int j=0;j<2*B_SZ*2*B_SZ;j++) {
				for (int i=0;i<2*B_SZ*2*B_SZ;i++) {
					covariance[j][i]/=scale[j]*scale[i];
				}
			}
		}

		protected int coeffDrop() {
			double delta=Double.MAX_VALUE;
			int y=-1;
			int c=-1;
			for (int j=0;j<3*B_SZ*B_SZ;j++) {
				for (int i=0;i<B_SZ*B_SZ;i++) {
					if (mult[j][i]) {
						if (mseMasked[j][i]-mse[i]<delta) {
							delta=mseMasked[j][i]-mse[i];
							y=i;
							c=j;
						}
					}
				}
			}
			mult[c][y]=false;
			return(y);
		}

		protected double mse(int _y) {
			return(mse(_y,false));
		}

		protected double mse(int _y,boolean _saveBeta) {
			double yty=covariance[3*B_SZ*B_SZ+_y][3*B_SZ*B_SZ+_y];
			double ytxb=0;
			int xsize=numMults[_y];
			if (xsize>0) {
				if (xtx[xsize-1]==null) {
					xty[xsize-1]=new DenseMatrix64F(xsize,1);
					b[xsize-1]=new DenseMatrix64F(xsize,1);
				}
				xtx[xsize-1]=new DenseMatrix64F(xsize,xsize);
				for (int ji=0,j=0;j<3*B_SZ*B_SZ;j++) {
					if (mult[j][_y]) {
						for (int ii=0,i=0;i<3*B_SZ*B_SZ;i++) {
							if (mult[i][_y]) {
								xtx[xsize-1].set(ji,ii,covariance[j][i]);
								ii++;
							}
						}
						xty[xsize-1].set(ji,0,covariance[j][3*B_SZ*B_SZ+_y]);
						ji++;
					}
				}

				solver.setA(xtx[xsize-1]);
				solver.solve(xty[xsize-1],b[xsize-1]);
				// compute y^T * x * beta
				for (int ii=0,i=0;i<3*B_SZ*B_SZ;i++) {
					if (mult[i][_y]) {
						ytxb+=covariance[3*B_SZ*B_SZ+_y][i]*b[xsize-1].get(ii,0);
						ii++;
					}
				}
				// save beta for use with classification
				if (_saveBeta) {
					for (int ii=0,i=0;i<3*B_SZ*B_SZ;i++) {
						if (mult[i][_y]) {
							// beta' = Sx^-1 * beta * Sy -> beta = Sx * beta' * Sy^-1 
							beta[i][_y]=b[xsize-1].get(ii,0)*scale[3*B_SZ*B_SZ+_y]/scale[i];
							ii++;
						}
						else {
							beta[i][_y]=0;
						}
					}
				}
			}
			return(yty-ytxb);
		}

		// compute the masked MSE's for a given y coefficient
		protected void updateMasked(int _y) {
			for (int i=0;i<3*B_SZ*B_SZ;i++) {
				if (mult[i][_y]) {
					mult[i][_y]=false;
					numMults[_y]--;
					mseMasked[i][_y]=mse(_y);
					mult[i][_y]=true;
					numMults[_y]++;
				 }
			 }
		}

		protected double predError(int[] _data,int _y) {
			double pred=0;
			for (int i=0;i<3*B_SZ*B_SZ;i++) {
				pred+=_data[INDEX[i]]*beta[i][_y];
			}
			return(Math.abs(_data[INDEX[3*B_SZ*B_SZ+_y]]-pred));
		}

		protected double mseUpdate(double _error,int _y,int _weight) {
			double se=_error*_error;
			double wt=weightTotal+_weight;
			double delta=se-mse[_y];
			return(mse[_y]+delta*_weight/wt);
		}

	};

	protected static LinearSolver<DenseMatrix64F> solver=new SolvePseudoInverseSvd();

	protected static DenseMatrix64F[] xtx=new DenseMatrix64F[3*B_SZ*B_SZ];

	protected static DenseMatrix64F[] xty=new DenseMatrix64F[3*B_SZ*B_SZ];

	protected static DenseMatrix64F[] b=new DenseMatrix64F[3*B_SZ*B_SZ];

	protected static double satdAverage(ModeData[] _modes) {
		// compute online update of mean
		return(0);
	}

	public static void main(String[] _args) throws IOException {
		ModeData[] modeData=new ModeData[MODES];

		// initialize the modes
		for (int i=0;i<MODES;i++) {
			modeData[i]=new ModeData();
		}

		// for each pass of k-means load the data
		File folder=new File("data");
		File[] coeffFiles=folder.listFiles(new FilenameFilter() {
			public boolean accept(File _file,String _filename) {
				return(_filename.endsWith(".coeffs"));
			}
		});
		if (coeffFiles==null) {
			System.out.println("No .coeffs files found in data/ folder.  Enable the PRINT_BLOCKS ifdef and run:");
			System.out.println("  for f in subset1-y4m/*.y4m; do ./init_intra_maps $f && ./init_intra_xform $f 2> $f.coeffs; done");
			System.out.println("Move the *.coeffs files to tools/java/data");
			System.exit(0);
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

		for (File file : binFiles) {
			System.out.println(file.getPath());
			DataInputStream dis=new DataInputStream(new BufferedInputStream(new FileInputStream(file),DISK_BLOCK_SIZE));
			DataOutputStream dos=new DataOutputStream(new BufferedOutputStream(new FileOutputStream(file.getPath()+".step00.mode"),DISK_BLOCK_SIZE));
			// 3 color planes, YUV
			for (int pli=0;pli<1;pli++) {
				// read the plane size
				int nx=dis.readShort();
				int ny=dis.readShort();
				System.out.println("nx="+nx+" ny="+ny);

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
						//modeData[mode].addBlock(mode==0?1:weight,data);
						modeData[mode].addBlock(weight,data);
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

		/*System.out.println("numBlocks="+modeData[0].numBlocks);
		for (int i=0;i<3*B_SZ*B_SZ;i++) {
			System.out.println(i+": "+modeData[0].mean[i]);
		}
		System.exit(0);*/

		long last=System.currentTimeMillis();
		for (int step=1;step<=30;step++) {
			System.out.println("*** Starting step "+step);

			// phase 1:
			//   normalize
			//   compute the metric (MSE or CG)
			//   sparsify (optional)
			//   compute beta

			long modeStart=System.currentTimeMillis();
			// for each mode
			for (int mode=0;mode<MODES;mode++) {
				ModeData md=modeData[mode];

				md.normalize();

				// compute MSE based on completely unmasked

				double mse=0;
				double cg=0;
				for (int y=0;y<B_SZ*B_SZ;y++) {
					md.mse[y]=md.mse(y,true);
					cg+=-10*Math.log10(md.mse[y]);
					mse+=md.mse[y];
				}

				System.out.println("  "+mode+": blocks="+md.numBlocks+"\tMSE="+mse+"\tCG="+cg);

				double metric=USE_CG?cg:mse;

				//System.out.println("before="+(System.currentTimeMillis()-modeStart));
				int mults=(3*B_SZ*B_SZ)*(B_SZ*B_SZ);
				if (DROP_MULTS) {
					for (int y=0;y<B_SZ*B_SZ;y++) {
						md.updateMasked(y);
					}
					//double rmseStart=Math.sqrt(md.total/(B_SZ*B_SZ));
					for (;mults-->MAX_MULTS;) {
						int y=modeData[mode].coeffDrop();
						metric-=USE_CG?-10*Math.log10(md.mse[y]):md.mse[y];
						md.mse[y]=md.mse(y,true);
						metric+=USE_CG?-10*Math.log10(md.mse[y]):md.mse[y];
						md.numMults[y]--;
						md.updateMasked(y);
					}
					//double rmseEnd=Math.sqrt(md.total/(B_SZ*B_SZ));
					/*System.out.println(mode+" "+rmseEnd+" "+(rmseEnd-rmseStart)+" took "+(modeEnd-modeStart)/1000.0);
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
					}*/
				}
				//System.out.println("after="+(System.currentTimeMillis()-modeStart));
				long modeEnd=System.currentTimeMillis();
				modeStart=modeEnd;
			}

			ModeData[] oldModeData=modeData;

			// reset the mode data
			modeData=new ModeData[MODES];
			for (int i=0;i<MODES;i++) {
				modeData[i]=new ModeData();
			}

			double total_satd=0;
			int total_blocks=0;

			// re-classify the blocks based on the computed beta's
			for (File file : binFiles) {
				//System.out.print(file.getPath());
				DataInputStream dis=new DataInputStream(new BufferedInputStream(new FileInputStream(file),DISK_BLOCK_SIZE));
				DataInputStream dis2=new DataInputStream(new BufferedInputStream(new FileInputStream(file.getPath()+".step"+(step-1<10?"0"+(step-1):step-1)+".mode"),DISK_BLOCK_SIZE));
				DataOutputStream dos=new DataOutputStream(new BufferedOutputStream(new FileOutputStream(file.getPath()+".step"+(step<10?"0"+step:step)+".mode"),DISK_BLOCK_SIZE));
				// 3 color planes, YUV
				for (int pli=0;pli<1;pli++) {
					// read the plane size
					int nx=dis.readShort();
					int ny=dis.readShort();

					int[] to=new int[MODES];
					int[] from=new int[MODES];

					double[] image_satd=new double[MODES];
					int[] image_blocks=new int[MODES];

					int[] rgb=new int[nx*ny];
					for (int block=0;block<nx*ny;block++) {
						// load the data
						int mode=dis.readShort();
						int weight=dis.readShort();
						int oldMode=dis2.readShort();
						int[] data=new int[2*B_SZ*2*B_SZ];
						for (int i=0;i<2*B_SZ*2*B_SZ;i++) {
							data[i]=dis.readShort();
						}
						//System.out.println("+++ block="+block+" mode="+mode+" weight="+weight);
						if (weight>0) {
							double[] errors=new double[MODES];
							double[] bits=new double[MODES];
							for (int j=0;j<MODES;j++) {
								ModeData md=oldModeData[j];
								for (int y=0;y<B_SZ*B_SZ;y++) {
									double error=md.predError(data,y);
									errors[j]+=error;
									//int tmpWeight=j==0?1:weight;
									bits[j]+=-10*Math.log10(md.mse[y]/md.mseUpdate(error,y,j==oldMode?-weight:weight));
									/*System.out.println("mode="+j+" y="+y);
									System.out.println("  old mse="+md.mse[y]);
									double oldBits=-10*Math.log10(md.mse[y]);
									System.out.println("  old bits="+oldBits);
									System.out.println("  error="+error);
									double mse=md.mseUpdate(error,y,-weight);
									System.out.println("  new mse="+mse);
									double newBits=-10*Math.log10(mse);
									System.out.println("  new bits="+newBits);
									double delta_bits=-10*Math.log10(md.mse[y]/mse);
									System.out.println("  delta bits="+delta_bits);*/
								}
							}
							double best=-Double.MAX_VALUE;
							double nextBest=-Double.MAX_VALUE;
							for (int j=0;j<MODES;j++) {
								double metric=0;
								if (USE_CG) {
									if (oldMode!=j) {
										// greater CG is better
										metric=bits[oldMode]+bits[j];
									}
								}
								else {
									// lower SATD is better
									metric=errors[oldMode]-errors[j];
								}
								/*for (int y=0;y<B_SZ*B_SZ;y++) {
									// compute |y-x*B| the SATD (or L1-norm)
									double error=md.predError(data,y);
									metric+=USE_CG?-10*Math.log10(md.mse[y]/md.mseUpdate(error,y,weight)):error;
								}
								//System.out.println("mode="+j+" metric="+metric);
								if (USE_CG) {
									metric=-(mode_bits+metric);
								}
								if (metric<best) {
									nextBest=best;
									best=metric;
									mode=j;
								}
								else {
									if (metric<nextBest) {
										nextBest=metric;
									}
								}*/
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
							//System.out.println("block="+block+" best="+best+" new_mode="+mode+" nextBest="+nextBest);
							image_satd[mode]+=errors[mode];
							image_blocks[mode]++;

							if (mode!=oldMode) {
								to[mode]++;
								from[oldMode]++;
							}
							modeData[mode].addBlock(mode==0?1:best-nextBest,data);
							//modeData[mode].addBlock(mode==0?1:weight,data);
							//modeData[mode].addBlock(weight,data);
						}
						dos.writeShort(mode);
						rgb[block]=MODE_COLORS[mode];
					}
					System.out.println("Image Average SATD");
					int total_to=0;
					int total_from=0;
					int blocks=0;
					for (int mode=0;mode<MODES;mode++) {
						System.out.println("  "+mode+": "+image_satd[mode]/image_blocks[mode]+"\t"+image_blocks[mode]+"\t-"+from[mode]+"\t+"+to[mode]);
						total_satd+=image_satd[mode];
						total_to+=to[mode];
						total_from+=from[mode];
						blocks+=image_blocks[mode];
					}
					System.out.println("\t\t\t"+blocks+"\t-"+total_from+"\t+"+total_to);
					//total_blocks+=blocks;
					total_blocks+=nx*ny;
					BufferedImage bi=new BufferedImage(nx,ny,BufferedImage.TYPE_INT_ARGB);
					bi.setRGB(0,0,nx,ny,rgb,0,nx);
					ImageIO.write(bi,"PNG",new File(file.getPath()+".step"+(step<10?"0"+step:step)+".png"));
				}
				dis.close();
				dis2.close();
				dos.close();
			}

			long now=System.currentTimeMillis();
			System.out.println("*** Ending step "+step+" took "+(now-last)/1000.0);
			//System.out.println("Step "+step+" took "+(now-last)/1000.0);
			System.out.println("Average SATD: "+(total_satd/total_blocks));
		}
	}

};
