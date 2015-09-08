/*
 * #%L
 * ImgLib2: a general-purpose, multidimensional image processing library.
 * %%
 * Copyright (C) 2009 - 2014 Stephan Preibisch, Tobias Pietzsch, Barry DeZonia,
 * Stephan Saalfeld, Albert Cardona, Curtis Rueden, Christian Dietz, Jean-Yves
 * Tinevez, Johannes Schindelin, Lee Kamentsky, Larry Lindsey, Grant Harris,
 * Mark Hiner, Aivar Grislis, Martin Horn, Nick Perry, Michael Zinsmaier,
 * Steffen Jaensch, Jan Funke, Mark Longair, and Dimiter Prodanov.
 * %%
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 * #L%
 */

package net.imglib2.algorithm.fft2;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;

import net.imglib2.Cursor;
import net.imglib2.Dimensions;
import net.imglib2.FinalDimensions;
import net.imglib2.Point;
import net.imglib2.algorithm.region.hypersphere.HyperSphere;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.type.numeric.complex.ComplexFloatType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.type.numeric.complex.ComplexFloatType;
import net.imglib2.type.numeric.real.FloatType;

import org.jtransforms.fft.FloatFFT_1D;

import edu.mines.jtk.dsp.FftComplex;
import edu.mines.jtk.dsp.FftReal;

import org.junit.Test;

import com.carrotsearch.junitbenchmarks.AbstractBenchmark;

/**
 * Tests {@link FFT}.
 * 
 * @author Curtis Rueden
 * @author Brian Northan
 */
public class FFTTest {

	/**
	 * Test basic FFT functionality by performing FFT followed by inverse FFT then
	 * asserting that the original image has been recovered
	 */
	@Test
	public void testFFT() {

		for (int i = 40; i < 90; i++) {
			long[] dim = new long[] { i, i, i };

			Img<FloatType> input = new ArrayImgFactory<FloatType>().create(dim,
				new FloatType());
			placeSphereInCenter(input);

			ImgFactory<ComplexFloatType> fftImgFactory = null;

			try {
				fftImgFactory = input.factory().imgFactory(new ComplexFloatType());
			}
			catch (IncompatibleTypeException ex) {
				fftImgFactory = null;
			}

			Img<ComplexFloatType> fft = FFT.realToComplex(input, fftImgFactory);

			Img<FloatType> inverse = input.factory().create(new FinalDimensions(dim),
				new FloatType());
			FFT.complexToRealUnpad(fft, inverse);

			assertImagesEqual(input, inverse, 0.00005f);
		}
	}

	/**
	 * confirm 2D jtransform is producing the same result as the mines version
	 */
	@Test
	public void testJTransformVsMines() {
		
		int trials=1000;

		for (int i = 50; i < 1048; i++) {

			Dimensions unpaddedDimensions = new FinalDimensions(new long[] { i, i });

			long[] paddedDimensionsMines = new long[2];
			long[] fftDimensionsMines = new long[2];

			// use the "mines" version of FFTMethods to calculate the dimensions
			net.imglib2.algorithm.fft2.mines.FFTMethods.dimensionsRealToComplexFast(
				unpaddedDimensions, paddedDimensionsMines, fftDimensionsMines);

			Dimensions dimensionsMines = new FinalDimensions(paddedDimensionsMines);

			long[] paddedDimensionsJTransform = new long[2];
			long[] fftDimensionsJTransform = new long[2];

			// use the 'jtransform' version of FFTMethods to confirm that
			// jtransform supports the dimensions
			FFTMethods.dimensionsRealToComplexSmall(dimensionsMines,
				paddedDimensionsJTransform, fftDimensionsJTransform);

			Dimensions dimensions = new FinalDimensions(paddedDimensionsMines);

			// assert that the mines dimensions and jtransform dimensions are
			// equal
			assertEquals(paddedDimensionsMines[0], paddedDimensionsJTransform[0]);
			assertEquals(paddedDimensionsMines[1], paddedDimensionsJTransform[1]);

			// using the above calculated dimensions to form an input array

			int size = (int) (paddedDimensionsMines[0] * paddedDimensionsMines[1]);
			final float[] array = new float[size];
			final float[] array2 = new float[size];

			// place a couple of impulses in the array
			array[size / 4] = 50;
			array[size / 2] = 100;

			// make a copy of the array
			for (int j = 0; j < size; j++) {
				array2[j] = array[j];
			}

			// create images for input and transform
			Img<FloatType> inJTransform = ArrayImgs.floats(array,
				paddedDimensionsMines);
			Img<FloatType> inMines = ArrayImgs.floats(array2, paddedDimensionsMines);

			Img<ComplexFloatType> transformJTransform =
				new ArrayImgFactory<ComplexFloatType>().create(fftDimensionsMines,
					new ComplexFloatType());
			Img<ComplexFloatType> transformMines =
				new ArrayImgFactory<ComplexFloatType>().create(fftDimensionsMines,
					new ComplexFloatType());
			
			FFTMethods.realToComplex(inJTransform, transformJTransform, 0);

			long start = System.currentTimeMillis();

			for (int j = 0; j < trials; j++) {
				// perform forward and inverse fft using the FFTMethods (which uses
				// jtransform)
				FFTMethods.realToComplex(inJTransform, transformJTransform, 0);
			}

			long jTransformTime = System.currentTimeMillis() - start;

			System.out.println();
			System.out.println("i: " + i);
			System.out.println("Jtransform time: " + jTransformTime);

			start = System.currentTimeMillis();

			for (int j = 0; j < trials; j++) {
				// perform forward and inverse fft with the legacy FFTMethods which
				// uses mines
				net.imglib2.algorithm.fft2.mines.FFTMethods.realToComplex(inMines,
					transformMines, 0);
			}

			long minesTime = System.currentTimeMillis() - start;

			System.out.println("Mines time: " + minesTime);
			System.out.println();

			// assert that the images are the same
			assertComplexImagesEqual(transformJTransform, transformMines, 0.0001f);

		}
	}

	@Test
	public void testFFTMethodsVsJTransformAPIBenchmark() {

		int startSize = 90;
		int finishSize = 110;
		int incr=1;
		
		// number of trials to perform for each size
		int trials = 1000;
		
		// a list of sizes for which jTransform is faster
		ArrayList<Integer> jTransformFaster=new ArrayList<Integer>();
		
		// a list of sizes for which jTransform is 10% of mines
		ArrayList<Integer> jTransformWithinX=new ArrayList<Integer>();
		Double X=1.10;
		
		for (int size = startSize; size < finishSize; size+=incr) {

			int jTransformSize = NextSmoothNumber.nextSmooth(size);
			float[] array = createTestArray(jTransformSize);

			int minesSize = FftReal.nfftFast(size);
			final int complexSize = (minesSize / 2 + 1) * 2;

			float[] minesArrayIn = createTestArray(minesSize);
			float[] minesArrayOut = new float[complexSize];

			// create and perform forward and inverse fft on the data using
			// JTransform
			// directly
			final FloatFFT_1D fft_jt = new FloatFFT_1D(jTransformSize);
			
			final FftReal fft_mines = new FftReal(minesSize);

			long start = System.currentTimeMillis();

			for (int i = 1; i < trials; i++) {
				minesArrayIn = createTestArray(minesSize);
				fft_mines.realToComplex(1, minesArrayIn, minesArrayOut);
				//fft_mines.complexToReal(1, minesArrayOut, minesArrayIn);
			}

			long minesTime = System.currentTimeMillis() - start;

			start = System.currentTimeMillis();
			
			fft_jt.realForward(array);
			fft_jt.realInverse(array, true);

			for (int i = 1; i < trials; i++) {
				array = createTestArray(jTransformSize);
				fft_jt.realForward(array);
				//fft_jt.realInverse(array, true);
			}

			long jTransformTime = System.currentTimeMillis() - start;

			System.out.println("size: " + size);
			System.out.println("jtransform size: "+jTransformSize+" mines size: "+minesSize);
			System.out.println("JTransform time: " + jTransformTime);
			System.out.println("Mines time: " + minesTime);
			
			if (jTransformTime<minesTime) {
				jTransformFaster.add(size);
			}
			else if (jTransformTime<X*minesTime) {
				jTransformWithinX.add(size);
			}
		}
		
		System.out.println("jTransform faster: ");
		System.out.println(jTransformFaster);
		
		System.out.println("jTransform within "+X);
		System.out.println(jTransformWithinX);
	}

	/**
	 * This test confirms that the result using FFTMethods is equal to the result
	 * obtained using the JTransform api
	 */
	@Test
	public void testFFTMethodsVsJTransformAPI() {

		for (int i = 40; i < 100; i++) {

			// determine the dimensions that FFTMethods will use given i
			Dimensions unpaddedDimensions = new FinalDimensions(new long[] { i });
			long[] paddedDimensions = new long[1];
			long[] fftDimensions = new long[1];

			FFTMethods.dimensionsRealToComplexSmall(unpaddedDimensions,
				paddedDimensions, fftDimensions);
			Dimensions dimensions = new FinalDimensions(paddedDimensions);

			// using the above calculated dimensions form an input array
			int size = (int) paddedDimensions[0];
			final float[] array = new float[size];
			final float[] copy = new float[size];

			// place a couple of impulses in the array
			array[size / 4] = 50;
			array[size / 2] = 100;

			// make a copy of the array
			for (int j = 0; j < size; j++) {
				copy[j] = array[j];
			}

			// create and perform forward and inverse fft on the data using
			// JTransform
			// directly
			final FloatFFT_1D fft = new FloatFFT_1D(size);

			fft.realForward(array);
			fft.realInverse(array, true);

			// assert that we get the original signal back (within an error
			// delta)
			for (int j = 0; j < size; j++) {
				assertEquals(copy[j], array[j], 0.001);
			}

			// now use the FFTMethods api

			// create images for input, transform and inverse
			Img<FloatType> in = ArrayImgs.floats(copy, paddedDimensions);

			Img<ComplexFloatType> transform = new ArrayImgFactory<ComplexFloatType>()
				.create(fftDimensions, new ComplexFloatType());

			Img<FloatType> inverse = new ArrayImgFactory<FloatType>().create(
				dimensions, new FloatType());

			// perform forward and inverse fft using the FFTMethods approach
			FFTMethods.realToComplex(in, transform, 0);
			FFTMethods.complexToReal(transform, inverse, 0);

			int j = 0;

			Cursor<FloatType> cin = in.cursor();
			Cursor<FloatType> cinverse = inverse.cursor();

			while (cin.hasNext()) {

				cin.fwd();
				cinverse.fwd();

				// assert that the inverse = the input within the error delta
				assertEquals(cin.get().getRealFloat(), cinverse.get().getRealFloat(),
					0.001);

				// assert that the inverse obtained using FFTMethods api is
				// exactly
				// equal to the inverse obtained from using JTransform directly
				assertEquals(array[j], cinverse.get().getRealFloat(), 0);

				j++;
			}

		}
	}

	/**
	 * Test 1D JTransform FFT and mines FFT
	 */
	@Test
	public void testJTransformVsMines1D() {

		int size = 11;

		// determine the dimensions that FFTMethods will use given i
		Dimensions unpaddedDimensions = new FinalDimensions(new long[] { size });

		long[] paddedDimensionsMines = new long[1];
		long[] fftDimensionsMines = new long[1];

		net.imglib2.algorithm.fft2.mines.FFTMethods.dimensionsRealToComplexFast(
			unpaddedDimensions, paddedDimensionsMines, fftDimensionsMines);

		// using the above calculated dimensions form an input array
		int realsize = (int) paddedDimensionsMines[0];
		int complexsize = (int) fftDimensionsMines[0];

		final float[] memJTransform = new float[realsize];
		final float[] memInMines = new float[realsize];
		final float[] memOutMines = new float[2 * complexsize];

		// place a couple of impulses in the array
		memJTransform[realsize / 4] = 50;
		memJTransform[realsize / 2] = 100;

		// make a copy of the array
		for (int j = 0; j < realsize; j++) {
			memInMines[j] = memJTransform[j];
		}

		// create and perform forward and inverse fft on the data using
		// JTransform
		// directly
		final FloatFFT_1D fft = new FloatFFT_1D(realsize);
		fft.realForward(memJTransform);

		// mines directly
		final FftReal fftmines = new FftReal(realsize);
		fftmines.realToComplex(-1, memInMines, memOutMines);

		// create images for input, transform
		Img<FloatType> inJTransform = ArrayImgs.floats(memInMines,
			paddedDimensionsMines);
		Img<ComplexFloatType> fftJTransform =
			new ArrayImgFactory<ComplexFloatType>().create(fftDimensionsMines,
				new ComplexFloatType());

		Img<FloatType> inMines = ArrayImgs.floats(memInMines,
			paddedDimensionsMines);
		Img<ComplexFloatType> fftMines = new ArrayImgFactory<ComplexFloatType>()
			.create(fftDimensionsMines, new ComplexFloatType());

		// perform forward and inverse fft using the FFTMethods approach
		FFTMethods.realToComplex(inJTransform, fftJTransform, 0);

		net.imglib2.algorithm.fft2.mines.FFTMethods.realToComplex(inMines, fftMines,
			0);

		for (int j = 0; j < realsize; j++) {
			System.out.println("j: " + j + " jt: " + memJTransform[j] + " mines: " +
				memOutMines[j]);

		}

		System.out.println("mt r: " + memOutMines[realsize]);
		System.out.println("mt c: " + memOutMines[realsize + 1]);

		Cursor<ComplexFloatType> cjt = fftJTransform.cursor();
		Cursor<ComplexFloatType> cmines = fftMines.cursor();

		// assert that the jtransform and mines results are the same
		while (cjt.hasNext()) {

			cjt.fwd();
			cmines.fwd();

			System.out.println();
			System.out.println("jt: " + cjt.get().getRealFloat() + " " + cjt.get()
				.getImaginaryFloat());
			System.out.println("mines: " + cmines.get().getRealFloat() + " " + cmines
				.get().getImaginaryFloat());

			assertEquals(cjt.get().getRealFloat(), cmines.get().getRealFloat(),
				0.0001);
			assertEquals(cjt.get().getImaginaryFloat(), cmines.get()
				.getImaginaryFloat(), 0.0001);
		}

		// finally also do a hard coded test

		cjt.reset();
		cjt.fwd();
		assertEquals(cjt.get().getRealFloat(), 150.0, 0.0001);
		assertEquals(cjt.get().getImaginaryFloat(), 0.0, 0.001);
		cjt.fwd();
		assertEquals(cjt.get().getRealFloat(), -100.0, 0.0001);
		assertEquals(cjt.get().getImaginaryFloat(), -50.0, 0.001);
		cjt.fwd();
		assertEquals(cjt.get().getRealFloat(), 50.0, 0.0001);
		assertEquals(cjt.get().getImaginaryFloat(), 3.2898595E-7, 0.001);
		cjt.fwd();
		assertEquals(cjt.get().getRealFloat(), -100.0, 0.0001);
		assertEquals(cjt.get().getImaginaryFloat(), 50.0, 0.001);
		cjt.fwd();
		assertEquals(cjt.get().getRealFloat(), 150.0, 0.0001);
		assertEquals(cjt.get().getImaginaryFloat(), 3.2898595E-7, 0.001);
		cjt.fwd();
		assertEquals(cjt.get().getRealFloat(), -100.0, 0.0001);
		assertEquals(cjt.get().getImaginaryFloat(), -50.0, 0.001);
		cjt.fwd();
		assertEquals(cjt.get().getRealFloat(), 50.0, 0.0001);
		assertEquals(cjt.get().getImaginaryFloat(), 0.0, 0.001);

	}

	/**
	 * a utility to assert that two images are equal
	 * 
	 * @param img1
	 * @param img2
	 * @param delta
	 */
	protected void assertImagesEqual(Img<FloatType> img1, Img<FloatType> img2,
		float delta)
	{
		Cursor<FloatType> c1 = img1.cursor();
		Cursor<FloatType> c2 = img2.cursor();

		while (c1.hasNext()) {

			c1.fwd();
			c2.fwd();

			// assert that the inverse = the input within the error delta
			assertEquals(c1.get().getRealFloat(), c2.get().getRealFloat(), delta);
		}

	}

	// a utility to assert that two images are equal
	protected void assertComplexImagesEqual(Img<ComplexFloatType> img1,
		Img<ComplexFloatType> img2, float delta)
	{
		Cursor<ComplexFloatType> c1 = img1.cursor();
		Cursor<ComplexFloatType> c2 = img2.cursor();

		int i = 0;
		while (c1.hasNext()) {

			c1.fwd();
			c2.fwd();

			i++;

			// assert that the inverse = the input within the error delta
			assertEquals(c1.get().getRealFloat(), c2.get().getRealFloat(), delta);
			// assert that the inverse = the input within the error delta
			assertEquals(c1.get().getImaginaryFloat(), c2.get().getImaginaryFloat(),
				delta);
		}

	}

	/**
	 * utility that places a sphere in the center of the image
	 * 
	 * @param img
	 */
	private void placeSphereInCenter(Img<FloatType> img) {

		final Point center = new Point(img.numDimensions());

		for (int d = 0; d < img.numDimensions(); d++)
			center.setPosition(img.dimension(d) / 2, d);

		HyperSphere<FloatType> hyperSphere = new HyperSphere<FloatType>(img, center,
			2);

		for (final FloatType value : hyperSphere) {
			value.setReal(1);
		}
	}

	private float[] createTestArray(int size) {
		final float[] array = new float[size];

		// place a couple of impulses in the array
		array[size / 4] = 50;
		array[size / 2] = 100;
		
		return array;
	}

}
