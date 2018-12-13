/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volvis;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;
import gui.RaycastRendererPanel;
import gui.TransferFunction2DEditor;
import gui.TransferFunctionEditor;
import java.awt.image.BufferedImage;
import util.TFChangeListener;
import util.VectorMath;
import volume.GradientVolume;
import volume.Volume;
import volume.VoxelGradient;

/**
 *
 * @author michel
 */
public class RaycastRenderer extends Renderer implements TFChangeListener {

    private Volume volume = null;
    private GradientVolume gradients = null;
    RaycastRendererPanel panel;
    TransferFunction tFunc;
    TransferFunctionEditor tfEditor;
    TransferFunction2DEditor tfEditor2D;

    //shading boolean
    public boolean shading = false;
    //rendering type variable
    public String TypeOfRender = "Slicer";

    public RaycastRenderer() {
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");
    }

    public void setVolume(Volume vol) {
        System.out.println("Assigning volume");
        volume = vol;

        System.out.println("Computing gradients");
        gradients = new GradientVolume(vol);

        // set up image for storing the resulting rendering
        // the image width and height are equal to the length of the volume diagonal
        int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX() * vol.getDimX() + vol.getDimY() * vol.getDimY()
                + vol.getDimZ() * vol.getDimZ()));
        if (imageSize % 2 != 0) {
            imageSize = imageSize + 1;
        }
        image = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_ARGB);
        // create a standard TF where lowest intensity maps to black, the highest to white, and opacity increases
        // linearly from 0.0 to 1.0 over the intensity range
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum());

        // uncomment this to initialize the TF with good starting values for the orange dataset 
        tFunc.setTestFunc();

        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram());

        tfEditor2D = new TransferFunction2DEditor(volume, gradients);
        tfEditor2D.addTFChangeListener(this);

        System.out.println("Finished initialization of RaycastRenderer");
    }

    public RaycastRendererPanel getPanel() {
        return panel;
    }

    public TransferFunction2DEditor getTF2DPanel() {
        return tfEditor2D;
    }

    public TransferFunctionEditor getTFPanel() {
        return tfEditor;
    }


    short getVoxel(double[] coord) {

        if (coord[0] < 0 || coord[0] >= volume.getDimX() || coord[1] < 0 || coord[1] >= volume.getDimY()
                || coord[2] < 0 || coord[2] >= volume.getDimZ()) {
            return 0;
        }

        //Interpolation
        //take floor values
        int x0 = (int) Math.floor(coord[0]);
        int y0 = (int) Math.floor(coord[1]);
        int z0 = (int) Math.floor(coord[2]);

        //take ceiling values
        int x1 = (int) Math.ceil(coord[0]);
        int y1 = (int) Math.ceil(coord[1]);
        int z1 = (int) Math.ceil(coord[2]);

        //calculate alpha, beta, gamma
        double alpha = 0;
        double beta = 0;
        double gamma = 0;
        if (coord[0] != x0) {
            alpha = (coord[0] - x0) / (x1 - x0);
        }

        if (coord[1] != y0) {
            beta = (coord[1] - y0) / (y1 - y0);
        }

        if (coord[2] != z0) {
            gamma = (coord[2] - z0) / (z1 - z0);
        }

        int val0 = volume.getVoxel(x0, y0, z0);
        int val1 = volume.getVoxel(x1, y0, z0);

        int val2 = volume.getVoxel(x0, y1, z0);
        int val3 = volume.getVoxel(x1, y1, z0);

        int val4 = volume.getVoxel(x0, y0, z1);
        int val5 = volume.getVoxel(x1, y0, z1);

        int val6 = volume.getVoxel(x0, y1, z1);
        int val7 = volume.getVoxel(x1, y1, z1);

        // interpolate along x
        double val01 = alpha * val1 + (1 - alpha) * val0;
        double val23 = alpha * val3 + (1 - alpha) * val2;
        double val45 = alpha * val5 + (1 - alpha) * val4;
        double val67 = alpha * val7 + (1 - alpha) * val6;

        // interpolate along y
        double val0123 = beta * val23 + (1 - beta) * val01;
        double val4567 = beta * val67 + (1 - beta) * val45;

        // interpolate along z and the final interpolated value
        double finalVal = gamma * val4567 + (1 - gamma) * val0123;

        return (short) Math.round(finalVal);
        //return volume.getVoxel(x, y, z);
    }

    VoxelGradient getVoxelGradient(double[] coord) {
        VoxelGradient voxelGradient = new VoxelGradient();

        int x = (int) Math.floor(coord[0]);
        int y = (int) Math.floor(coord[1]);
        int z = (int) Math.floor(coord[2]);

        if (coord[0] < 0 || coord[0] >= volume.getDimX() - 1 || coord[1] < 0 || coord[1] >= volume.getDimY() - 1
                || coord[2] < 0 || coord[2] >= volume.getDimZ() - 1) {
            return voxelGradient;
        }
        //Interpolation
        //take floor values
        int x0 = (int) Math.floor(coord[0]);
        int y0 = (int) Math.floor(coord[1]);
        int z0 = (int) Math.floor(coord[2]);

        //take ceiling values
        int x1 = (int) Math.ceil(coord[0]);
        int y1 = (int) Math.ceil(coord[1]);
        int z1 = (int) Math.ceil(coord[2]);


        //calculate alpha, beta, gamma
        float alpha = 0;
        float beta = 0;
        float gamma = 0;

        if (coord[0] != x0) {
            alpha = (float) (coord[0] - x0) / (x1 - x0);
        }

        if (coord[1] != y0) {
            beta = (float) (coord[1] - y0) / (y1 - y0);
        }

        if (coord[2] != z0) {
            gamma = (float) (coord[2] - z0) / (z1 - z0);
        }

        VoxelGradient val0 = gradients.getGradient(x0, y0, z0);
        VoxelGradient val1 = gradients.getGradient(x1, y0, z0);

        VoxelGradient val2 = gradients.getGradient(x0, y1, z0);
        VoxelGradient val3 = gradients.getGradient(x1, y1, z0);

        VoxelGradient val4 = gradients.getGradient(x0, y0, z1);
        VoxelGradient val5 = gradients.getGradient(x1, y0, z1);

        VoxelGradient val6 = gradients.getGradient(x0, y1, z1);
        VoxelGradient val7 = gradients.getGradient(x1, y1, z1);

        // interpolate along x
        VoxelGradient val01 = Add(Mult(alpha, val1), Mult((1 - alpha), val0));
        VoxelGradient val23 = Add(Mult(alpha, val3), Mult((1 - alpha), val2));
        VoxelGradient val45 = Add(Mult(alpha, val5), Mult((1 - alpha), val4));
        VoxelGradient val67 = Add(Mult(alpha, val7), Mult((1 - alpha), val6));

        // interpolate along y
        VoxelGradient val0123 = Add(Mult(beta, val23), Mult((1 - beta), val01));
        VoxelGradient val4567 = Add(Mult(beta, val67), Mult((1 - beta), val45));

        // interpolate along z and the final interpolated value
        VoxelGradient finalVal = Add(Mult(gamma, val4567), Mult((1 - gamma), val0123));

        return finalVal;

    }

    //multiply a 3-D voxelgradient with scalar value alpha
    public VoxelGradient Mult(float alpha, VoxelGradient val) {
        return new VoxelGradient(alpha * val.x, alpha * val.y, alpha * val.z);
    }

    //add two 3-D voxelgradient
    public VoxelGradient Add(VoxelGradient val1, VoxelGradient val2) {
        return new VoxelGradient(val1.x + val2.x, val1.y + val2.y, val1.z + val2.z);
    }

    //Apply Levoy method in calculating opacity
    public TFColor CalcLevoy(double[] voxelCoord, VoxelGradient grad, TFColor voxelColor) {
        TransferFunction2DEditor.TriangleWidget t = tfEditor2D.triangleWidget;
        int voxelIntensity = getVoxel(voxelCoord);
        voxelColor.r = t.color.r;
        voxelColor.g = t.color.g;
        voxelColor.b = t.color.b;

        // opacity reweighting
        if (grad.mag == 0 && voxelIntensity == t.baseIntensity) {
            voxelColor.a = t.color.a;
        }
        else if (grad.mag > 0 && ((voxelIntensity - t.radius * grad.mag) <= t.baseIntensity)
                && ((voxelIntensity + t.radius * grad.mag) >= t.baseIntensity)) {
            voxelColor.a = t.color.a * (1.0 - 1.0 / t.radius * (Math.abs(voxelIntensity - t.baseIntensity) / (grad.mag)));
        } else {
            voxelColor.a = 0;
        }
        return voxelColor;
    }

    //Apply Phong shading method
    public TFColor CalcPhong(VoxelGradient grad, double[] viewVec, TFColor voxelColor){
        double [] N = new double[3];
        double[] L = new double[3];
        double[] H = new double[3];
        double[] N1 = new double[3];
        double dotProductNL = 0;
        double dotProductVR = 0;

        //calculate N
        N[0] = grad.x/grad.mag;
        N[1] = grad.y/grad.mag;
        N[2] = grad.z/grad.mag;

        //calculate L
        VectorMath.setVector(L, -viewVec[0], -viewVec[1], -viewVec[2]);

        //calculate dot product of N and L
        dotProductNL = VectorMath.dotproduct(N, L);

         //calculate reflection vector R
        for (int i = 0; i < 3; i++) {
            N1[i] = N[i] *2 ;
        }
        double R1 = VectorMath.dotproduct(N1, L);
        double[] R = new double[3];

        for (int i = 0; i < 3; i++) {
            R[i] = R1 * N[i];
        }
        for (int i = 0; i < 3; i++) {
            R[i] = R[i] - L[i];
        }

        //calculate dot product between V and R
        dotProductVR = VectorMath.dotproduct(viewVec, R);

        //ca;culate the final values
        if (dotProductNL > 0 && dotProductVR > 0) {
            voxelColor.r = (0.1 * voxelColor.r) + (dotProductNL * 0.7 * voxelColor.r) + (voxelColor.r * 0.2 * (Math.pow(dotProductVR, 10)));
            voxelColor.g = (0.1 * voxelColor.g) + (dotProductNL * 0.7 * voxelColor.g) + (0.2 * voxelColor.g * (Math.pow(dotProductVR, 10)));
            voxelColor.b = (0.1 * voxelColor.b) + (dotProductNL * 0.7 * voxelColor.b) + (0.2 * voxelColor.b * (Math.pow(dotProductVR, 10)));
        }

        return voxelColor;
    }

    //Maximum Intesnity Projection
    void mip(double[] viewMatrix) {
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }
        // vector uVec and vVec define a plane through the origin,
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] voxelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);


        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();

        // interactive mode
        int sampleStep = 1;
        if (interactiveMode == true) {
            sampleStep = 5;
        }

        // define variables to control sample step length
        double XStep = viewVec[0] * sampleStep;
        double YStep = viewVec[1] * sampleStep;
        double ZStep = viewVec[2] * sampleStep;

        TFColor pixelColor = new TFColor();

        for (int j = 0; j < image.getHeight(); j++) {

            double voxelCoordXStart = uVec[0] * (-imageCenter) + vVec[0] * (j - imageCenter) + volumeCenter[0];
            double voxelCoordYStart = uVec[1] * (-imageCenter) + vVec[1] * (j - imageCenter) + volumeCenter[1];
            double voxelCoordZStart = uVec[2] * (-imageCenter) + vVec[2] * (j - imageCenter) + volumeCenter[2];

            for (int i = 0; i < image.getWidth(); i++) {

                boolean isIntersected = true;

                if(voxelCoordXStart < 0 || voxelCoordXStart >= volume.getDimX()) {
                        isIntersected = false;
                    }

                if(voxelCoordYStart < 0 || voxelCoordYStart >= volume.getDimY()) {
                        isIntersected = false;
                    }

                if(voxelCoordZStart < 0 || voxelCoordZStart >= volume.getDimZ()) {
                        isIntersected = false;
                    }

                int maxVoxelIntensity = 0;

                //take start and end
                long start = Math.round(-voxelCoordZStart);
                long end = Math.round((volume.getDimZ() - voxelCoordZStart));

                // ray intersects with volume
                if(isIntersected && start < end) {
                    voxelCoord[0] = voxelCoordXStart + end * viewVec[0];
                    voxelCoord[1] = voxelCoordYStart + end * viewVec[1];
                    voxelCoord[2] = voxelCoordZStart + end * viewVec[2];

                    for(long u = start; u < end; u += sampleStep){
                        int voxelIntensity = getVoxel(voxelCoord);
                        if(voxelIntensity > maxVoxelIntensity){
                            maxVoxelIntensity = voxelIntensity;
                        }
                        // move to next sample voxel
                        voxelCoord[0] -= XStep;
                        voxelCoord[1] -= YStep;
                        voxelCoord[2] -= ZStep;
                    }
                    // Map the intensity to a grey value by linear scaling
                    pixelColor.r = maxVoxelIntensity / max;
                    pixelColor.g = pixelColor.r;
                    pixelColor.b = pixelColor.r;
                    pixelColor.a = maxVoxelIntensity > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque

                    // BufferedImage expects a pixel color packed as ARGB in an int

                    int c_alpha = pixelColor.a <= 1.0 ? (int) Math.floor(pixelColor.a * 255) : 255;
                    int c_red = pixelColor.r <= 1.0 ? (int) Math.floor(pixelColor.r * 255) : 255;
                    int c_green = pixelColor.g <= 1.0 ? (int) Math.floor(pixelColor.g * 255) : 255;
                    int c_blue = pixelColor.b <= 1.0 ? (int) Math.floor(pixelColor.b * 255) : 255;
                    int finalPixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                    image.setRGB(i, j, finalPixelColor);
                }
                else {
                    image.setRGB(i, j, 1 << 24);
                }

                // move to the next pixel
                voxelCoordXStart += uVec[0];
                voxelCoordYStart += uVec[1];
                voxelCoordZStart += uVec[2];

            }
        }
    }

    //Apply ray compositing
    void compositing(double[] viewMatrix) {
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                //image.setRGB(i, j, 0);
            }
        }
        // vector uVec and vVec define a plane through the origin,
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] voxelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // interactive mode
        int sampleStep = 1;
        if(interactiveMode == true) {
            sampleStep = 5;
        }

        double XStep = viewVec[0] * sampleStep;
        double YStep = viewVec[1] * sampleStep;
        double ZStep = viewVec[2] * sampleStep;

        TFColor tfColor = null;
        TFColor voxelColor = new TFColor();

        for (int j = 0; j < image.getHeight(); j ++) {

            double voxelCoordXStart = uVec[0] * (-imageCenter) + vVec[0] * (j - imageCenter) + volumeCenter[0];
            double voxelCoordYStart = uVec[1] * (-imageCenter) + vVec[1] * (j - imageCenter) + volumeCenter[1];
            double voxelCoordZStart = uVec[2] * (-imageCenter) + vVec[2] * (j - imageCenter) + volumeCenter[2];

            for (int i = 0; i < image.getWidth(); i++) {

                TFColor pixelColor = new TFColor(0, 0, 0, 1);

                // compute intersections
                boolean isIntersected = true;

                if(voxelCoordXStart < 0 || voxelCoordXStart >= volume.getDimX()) {
                        isIntersected = false;
                    }

                if(voxelCoordYStart < 0 || voxelCoordYStart >= volume.getDimY()) {
                        isIntersected = false;
                    }

                if(voxelCoordZStart < 0 || voxelCoordZStart >= volume.getDimZ()) {
                        isIntersected = false;
                    }

                //take start and end
                long start = Math.round(-voxelCoordZStart);;
                long end = Math.round((volume.getDimZ() - voxelCoordZStart));

                // ray intersects with volume
                if(isIntersected && start < end) {
                    voxelCoord[0] = voxelCoordXStart + end * viewVec[0];
                    voxelCoord[1] = voxelCoordYStart + end * viewVec[1];
                    voxelCoord[2] = voxelCoordZStart + end * viewVec[2];

                    for(long u = start; u < end; u += sampleStep){

                        // use Tansfer function to tansfer an intensity to a color
                        tfColor = tFunc.getColor(getVoxel(voxelCoord));
                        voxelColor.r = tfColor.r;
                        voxelColor.g = tfColor.g;
                        voxelColor.b = tfColor.b;
                        voxelColor.a = tfColor.a;

                        if(shading) {
                            voxelColor = CalcPhong(getVoxelGradient(voxelCoord), viewVec, voxelColor);
                        }
                        pixelColor.r = (1 - voxelColor.a) * pixelColor.r + voxelColor.a * voxelColor.r;
                        pixelColor.g = (1 - voxelColor.a) * pixelColor.g + voxelColor.a * voxelColor.g;
                        pixelColor.b = (1 - voxelColor.a) * pixelColor.b + voxelColor.a * voxelColor.b;
                        pixelColor.a = (1 - voxelColor.a) * pixelColor.a + voxelColor.a;

                        // move to the next sample
                        voxelCoord[0] -= XStep;
                        voxelCoord[1] -= YStep;
                        voxelCoord[2] -= ZStep;

                    }
                    // BufferedImage expects a pixel color packed as ARGB in an int
                    int c_alpha = pixelColor.a <= 1.0 ? (int) Math.floor(pixelColor.a * 255) : 255;
                    int c_red = pixelColor.r <= 1.0 ? (int) Math.floor(pixelColor.r * 255) : 255;
                    int c_green = pixelColor.g <= 1.0 ? (int) Math.floor(pixelColor.g * 255) : 255;
                    int c_blue = pixelColor.b <= 1.0 ? (int) Math.floor(pixelColor.b * 255) : 255;
                    int finalPixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                    image.setRGB(i, j, finalPixelColor);
                }
                else {
                    image.setRGB(i, j, 1 << 24);
                }

                // move to the next pixel
                voxelCoordXStart += uVec[0];
                voxelCoordYStart += uVec[1];
                voxelCoordZStart += uVec[2];
            }
        }
    }

    //Apply 2D transfer function
    void twoDimTransfer(double[] viewMatrix) {
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin,
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] voxelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX()/2, volume.getDimY()/2, volume.getDimZ()/2);

        // sample on a plane through the origin of the volume data
        TFColor voxelColor = new TFColor();
        short voxelIntensity = 0;
        VoxelGradient voxelGradient = new VoxelGradient();

        // interactive mode
        int sampleStep = 1;
        if(interactiveMode == true) {
            sampleStep = 5;
        }

        final double XStep = viewVec[0] * sampleStep;
        final double YStep = viewVec[1] * sampleStep;
        final double ZStep = viewVec[2] * sampleStep;

        for (int j = 0; j < image.getHeight(); j ++) {

            double voxelCoordXStart = uVec[0] * (-imageCenter) + vVec[0] * (j - imageCenter) + volumeCenter[0];
            double voxelCoordYStart = uVec[1] * (-imageCenter) + vVec[1] * (j - imageCenter) + volumeCenter[1];
            double voxelCoordZStart = uVec[2] * (-imageCenter) + vVec[2] * (j - imageCenter) + volumeCenter[2];

            for (int i = 0; i < image.getWidth(); i++) {

                TFColor pixelColor = new TFColor(0, 0, 0, 1);

                // compute intersections
                boolean isIntersected = true;
                if(voxelCoordXStart < 0 || voxelCoordXStart >= volume.getDimX()) {
                        isIntersected = false;
                }

                if(voxelCoordYStart < 0 || voxelCoordYStart >= volume.getDimY()) {
                        isIntersected = false;
                }


                if(voxelCoordZStart < 0 || voxelCoordZStart >= volume.getDimZ()) {
                        isIntersected = false;
                }

                long start = Math.round(-voxelCoordZStart);
                long end = Math.round(volume.getDimZ() - voxelCoordZStart);

                if(start > end) {
                    isIntersected = false;
                }

                // ray intersects with volume
                if(isIntersected) {
                    voxelCoord[0] = voxelCoordXStart + (end + sampleStep) * viewVec[0];
                    voxelCoord[1] = voxelCoordYStart + (end + sampleStep) * viewVec[1];
                    voxelCoord[2] = voxelCoordZStart + (end + sampleStep) * viewVec[2];

                    for(long u = start; u < end; u += sampleStep) {

                        // move to the next sample
                        voxelCoord[0] -= XStep;
                        voxelCoord[1] -= YStep;
                        voxelCoord[2] -= ZStep;

                        //implement Levoy's method
                        voxelColor = CalcLevoy(voxelCoord, getVoxelGradient(voxelCoord), voxelColor);

                        //implement Phong method
                        if(shading) {
                            voxelColor = CalcPhong(getVoxelGradient(voxelCoord), viewVec, voxelColor);
                        }

                        // composite voxel colors
                        pixelColor.r = (1 - voxelColor.a) * pixelColor.r + voxelColor.a * voxelColor.r;
                        pixelColor.g = (1 - voxelColor.a) * pixelColor.g + voxelColor.a * voxelColor.g;
                        pixelColor.b = (1 - voxelColor.a) * pixelColor.b + voxelColor.a * voxelColor.b;
                        pixelColor.a = (1 - voxelColor.a) * pixelColor.a + voxelColor.a;
                    }

                    // BufferedImage expects a pixel color packed as ARGB in an int
                    int c_alpha = pixelColor.a <= 1.0 ? (int) Math.floor(pixelColor.a * 255) : 255;
                    int c_red = pixelColor.r <= 1.0 ? (int) Math.floor(pixelColor.r * 255) : 255;
                    int c_green = pixelColor.g <= 1.0 ? (int) Math.floor(pixelColor.g * 255) : 255;
                    int c_blue = pixelColor.b <= 1.0 ? (int) Math.floor(pixelColor.b * 255) : 255;
                    int finalPixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                    image.setRGB(i, j, finalPixelColor);
                }
                else {
                    image.setRGB(i, j, 1 << 24);
                }

                // move to the next pixel
                voxelCoordXStart += uVec[0];
                voxelCoordYStart += uVec[1];
                voxelCoordZStart += uVec[2];
            }
        }
    }

        //Existing slicer method
        void slicer(double[] viewMatrix) {
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }
        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);


        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        System.out.println(volume.getDimX() / 2);
        System.out.println(volume.getDimY() / 2);
        System.out.println(volume.getDimZ() / 2);
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();

        
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2];


                int val = getVoxel(pixelCoord);
                
                // Map the intensity to a grey value by linear scaling
                voxelColor.r = val/max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        }
    }


    private void drawBoundingBox(GL2 gl) {
        gl.glPushAttrib(GL2.GL_CURRENT_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glColor4d(1.0, 1.0, 1.0, 1.0);
        gl.glLineWidth(1.5f);
        gl.glEnable(GL.GL_LINE_SMOOTH);
        gl.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glDisable(GL.GL_LINE_SMOOTH);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL2.GL_LIGHTING);
        gl.glPopAttrib();

    }

    @Override
    public void visualize(GL2 gl) {


        if (volume == null) {
            return;
        }

        drawBoundingBox(gl);
        gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);

        long startTime = System.currentTimeMillis();

        if ("Slicer".equals(TypeOfRender)){
            slicer(viewMatrix);
        }

        if ("MIP".equals(TypeOfRender)){
            mip(viewMatrix);
        }

        if ("Composite".equals(TypeOfRender)){
            compositing(viewMatrix);
        }

        if("2DTransfer".equals(TypeOfRender)){
            twoDimTransfer(viewMatrix);
        }
        
        long endTime = System.currentTimeMillis();
        double runningTime = (endTime - startTime);
        panel.setSpeedLabel(Double.toString(runningTime));

        Texture texture = AWTTextureIO.newTexture(gl.getGLProfile(), image, false);

        gl.glPushAttrib(GL2.GL_LIGHTING_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        // draw rendered image as a billboard texture
        texture.enable(gl);
        texture.bind(gl);
        double halfWidth = image.getWidth() / 2.0;
        gl.glPushMatrix();
        gl.glLoadIdentity();
        gl.glBegin(GL2.GL_QUADS);
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glTexCoord2d(0.0, 0.0);
        gl.glVertex3d(-halfWidth, -halfWidth, 0.0);
        gl.glTexCoord2d(0.0, 1.0);
        gl.glVertex3d(-halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 1.0);
        gl.glVertex3d(halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 0.0);
        gl.glVertex3d(halfWidth, -halfWidth, 0.0);
        gl.glEnd();
        texture.disable(gl);
        texture.destroy(gl);
        gl.glPopMatrix();

        gl.glPopAttrib();


        if (gl.glGetError() > 0) {
            System.out.println("some OpenGL error: " + gl.glGetError());
        }

    }
    private BufferedImage image;
    private double[] viewMatrix = new double[4 * 4];

    @Override
    public void changed() {
        for (int i=0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }
}
