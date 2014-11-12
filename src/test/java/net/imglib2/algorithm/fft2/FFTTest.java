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
import net.imglib2.Cursor;
import net.imglib2.Dimensions;
import net.imglib2.FinalDimensions;
import net.imglib2.Point;
import net.imglib2.algorithm.region.hypersphere.HyperSphere;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.type.numeric.complex.ComplexFloatType;
import net.imglib2.type.numeric.real.FloatType;

import org.jtransforms.fft.FloatFFT_1D;

import org.junit.Test;

import com.carrotsearch.junitbenchmarks.AbstractBenchmark;

/**
 * Tests {@link FFT}.
 * 
 * @author Curtis Rueden
 * @author Brian Northan
 */
public class FFTTest extends AbstractBenchmark {

	/**
	 * Test basic FFT functionality by performing FFT followed by inverse FFT then
	 * asserting that the original image has been recovered
	 */
	@Test
	public void testFFT() {

		for (int i = 40; i < 90; i++) {
			long[] dim = new long[] { i, i, i };

			Img<FloatType> input =
				new ArrayImgFactory<FloatType>().create(dim, new FloatType());
			placeSphereInCenter(input);

			ImgFactory<ComplexFloatType> fftImgFactory = null;

			try {
				fftImgFactory = input.factory().imgFactory(new ComplexFloatType());
			}
			catch (IncompatibleTypeException ex) {
				fftImgFactory = null;
			}

			Img<ComplexFloatType> fft = FFT.realToComplex(input, fftImgFactory);

			Img<FloatType> inverse =
				input.factory().create(new FinalDimensions(dim), new FloatType());
			FFT.complexToRealUnpad(fft, inverse);

			assertImagesEqual(input, inverse, 0.00001f);
		}
	}

	/**
	 * Test 1D FFT using JTransform directly and using FFTMethods interface
	 */
	@Test
	public void testJTransformFFT1D() {

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

			// create and perform forward and inverse fft on the data using JTransform
			// directly
			final FloatFFT_1D fft = new FloatFFT_1D(size);

			fft.realForward(array);
			fft.realInverse(array, true);

			// assert that we get the original signal back (within an error delta)
			for (int j = 0; j < size; j++) {
				assertEquals(copy[j], array[j], 0.001);
			}

			// now use the FFTMethods api

			// create images for input, transform and inverse
			Img<FloatType> in = ArrayImgs.floats(copy, paddedDimensions);

			Img<ComplexFloatType> transform =
				new ArrayImgFactory<ComplexFloatType>().create(fftDimensions,
					new ComplexFloatType());

			Img<FloatType> inverse =
				new ArrayImgFactory<FloatType>().create(dimensions, new FloatType());

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

				// assert that the inverse obtained using FFTMethods api is exactly
				// equal to the inverse obtained from using JTransform directly
				assertEquals(array[j], cinverse.get().getRealFloat(), 0);
				j++;
			}
		}
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

		int i = 0;
		while (c1.hasNext()) {

			c1.fwd();
			c2.fwd();

			i++;

			// assert that the inverse = the input within the error delta
			assertEquals(c1.get().getRealFloat(), c2.get().getRealFloat(), delta);
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

		HyperSphere<FloatType> hyperSphere =
			new HyperSphere<FloatType>(img, center, 2);

		for (final FloatType value : hyperSphere) {
			value.setReal(1);
		}
	}

}
