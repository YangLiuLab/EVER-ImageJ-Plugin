import ij.plugin.*;
import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.filter.*;
import ij.io.Opener;
import java.io.*;

public class EVERv3_ implements PlugIn {

	private static float Baseline=(float)100, Gain=(float)0.47;
	private static int BGkernelwidth=(int)0;
	private static float BGoffset=(float)100;
	private static boolean flag=false;

	public void run(String arg) {
		ImagePlus imp = IJ.getImage();
		InputStream is = getClass().getResourceAsStream("EVER_LUT.tif");
		Opener opener = new Opener();
		ImagePlus implut =opener.openTiff(is, "EVER_LUT");
		
 		if (!showDialog())
			return;
		BGestimate(imp, Baseline, Gain, BGkernelwidth, BGoffset, flag, implut);

	}
	
	private boolean showDialog() {
		GenericDialog gd = new GenericDialog("EVER background correction");
		gd.addNumericField("Camera baseline (counts) : ", Baseline, 2);
		gd.addNumericField("Camera gain (photons/count) : ", Gain, 2);
		gd.addNumericField("Background smooth kernel width (pixel) : ", BGkernelwidth, 0);
		gd.addNumericField("Background subtraction offset (photons): ", BGoffset, 0);
		gd.addCheckbox("Generate background Stack", flag);
		gd.showDialog();
		if (gd.wasCanceled())
			return false;
		Baseline= (float) gd.getNextNumber();
		Gain=  (float) gd.getNextNumber();
		BGkernelwidth = (int) gd.getNextNumber();
		BGoffset = (float) gd.getNextNumber();
		flag= (boolean) gd.getNextBoolean();
		return true;
	}
	
	private void BGestimate(ImagePlus imp, float Baseline, float Gain, int BGkernelwidth, float BGoffset, boolean flag, ImagePlus implut) {
		ImageStack stack1 = imp.getImageStack();
		ImageStack stack = stack1.convertToFloat();
		ImageStack newstack = new ImageStack(stack.getWidth(), stack.getHeight());
		int imW = stack.getWidth();
		int imH = stack.getHeight();
		int dimension = imW*imH ;
		
		float[] pixels = new float[dimension];
		float[] imMIN = new float[dimension];
		float[] imMEN = new float[dimension];
		float[] imSTD = new float[dimension];
		float[] temp = new float[dimension];
		float[] imAVG = new float[100];
		float[] imRatio = new float[100];
		float imRatioAVG = (float)0;
		
		ImageStack bglut = implut.getImageStack();
		float[] BG_LUT = new float[bglut.getWidth()*bglut.getHeight()];
		BG_LUT = (float[]) bglut.getPixels(1);
		int lutW = bglut.getWidth();
		
		
		for (int f=1;f<=stack.getSize();f++) {
			pixels = (float[]) stack.getPixels(f);
			for (int j=0;j<dimension;j++) {
				pixels[j] = (float)((pixels[j] - Baseline)*Gain);
			}
		}
		
		for (int f=1;f<=100;f++) {
			pixels = (float[]) stack.getPixels(f);
			for (int j=0;j<dimension;j++) {
				imMEN[j] += (float)(pixels[j]);
				if((imMIN[j]>pixels [j])|(imMIN[j]<1))
					imMIN[j] = (float)pixels [j];
			}
		}

		for (int f=1;f<=100;f++) {
			pixels = (float[]) stack.getPixels(f);
			for (int j=0;j<dimension;j++) {
				imSTD[j] += (float)((pixels[j] - imMEN[j]/100)*(pixels[j] - imMEN[j]/100));
			}
		}
		
		double imBG0 =0;
		for (int j=0;j<dimension;j++) {
			imBG0 = 8.2469e-08*imMIN[j]*imMIN[j]*imMIN[j]-1.6481e-04*imMIN[j]*imMIN[j]+1.1546*imMIN[j]+13.3655;
			if(Math.sqrt(imSTD[j]/100)>2*Math.sqrt(imBG0))
				imSTD[j] =  0;
			else
				imSTD[j] = 1;
		}
		
		for(int fs=0;fs<stack.getSize()/100;fs++){
			imMIN = new float[dimension];
			temp = new float[dimension];
			imAVG = new float[100];
			imRatio = new float[100];
			imRatioAVG  =  (float)0;
			
			pixels = (float[]) stack.getPixels(1+fs*100);
			for (int j=0;j<dimension;j++) {
				imAVG[0] += (float)(pixels[j]*imSTD[j]);
				imMIN[j] = (float)pixels[j];
			}
			
			for (int f=2;f<=100;f++) {
				pixels = (float[])stack.getPixels(f+fs*100);
				for (int j=0;j<dimension;j++) {
					imAVG[f-1] += (float)(pixels[j]*imSTD[j]);
					if(imMIN[j]>pixels[j])
						imMIN[j] = (float)pixels[j];
				}

			}
			
			float imAVGmin = (float)imAVG[0];
			for (int n=0;n<100;n++) {
				imAVGmin = (imAVG[n]>imAVGmin)? imAVGmin : imAVG[n];
			}
			for (int n=0;n<100;n++) {
				imRatio[n] = imAVG[n]/imAVGmin ;
			}
			for (int n=0;n<100;n++) {
				imRatioAVG += imRatio[n]/100;
			}
			imRatioAVG = (imRatioAVG-1)*2 + 1;
			
			for (int i=0;i<imH;i++) 
				for (int j=0;j<imW;j++) {
					float ksum = (float)0;
					float knum = (float)0;
					for (int m=-BGkernelwidth;m<=BGkernelwidth;m++) 
						for (int n=-BGkernelwidth;n<=BGkernelwidth;n++) {
							if (((i+m)>-1) & ((i+m)<imH) & ((j+n)>-1) & ((j+n)<imW)) {
								knum  += 1.0;
								ksum  += (float) (imMIN[(i+m)*imW+j+n]);
							}
					}
					temp[i*imW+j] = ksum/knum;
			}
			for (int i=0;i<imH;i++) 
				for (int j=0;j<imW;j++) {
					imMIN[i*imW+j] = (float)temp[i*imW+j];
			}

			float T = 0;
			for (int f=1;f<=100;f++) {
				pixels = (float[]) stack.getPixels(f+fs*100);
				for (int j=0;j<dimension;j++) {
					T = (float)((BG_LUT[Math.round((imRatioAVG-1)*2000) * lutW + Math.round(imMIN[j])])*imRatio[f-1]);
					if (flag)
						pixels [j] = T;
					else if (pixels[j]>(T-BGoffset))
						pixels [j] = (float)(pixels[j]-(T-BGoffset)) ;
					else
						pixels [j] = (float)0;
				}
				newstack .addSlice(""+(f+fs*100), pixels );
			}
		}
		
		String stackname;
		if (flag)
			stackname = "Estimated Background";
		else
			stackname = "Background Corrected Image";
		ImagePlus nimp = new ImagePlus(stackname, newstack );
		imp.close();
		nimp.show();
	}
}
