package audio;

 /**
 * AudioManipulation.java
 *
 * Time-stamp: <2019-02-12 08:57:09 rlc3>
 *
 * Defines mixer/effect functions on audio streams
 * Utilises the AudioInputStream class 
 * 
 * To compile: javac -classpath editor.jar:. RunEffects.java
 * To run use: java -classpath editor.jar:. RunEffects
 * 
 */ 

import javax.sound.sampled.*;
import java.io.*;

public class AudioManipulation {

	/**** echo *****************************************************************/

	public static AudioInputStream echo(AudioInputStream ais, int timeDelay, double fading0, double fading1) {

		byte[] a = null;
		int[] data, ch0, ch1;
		int max;

		try {

			// AudioInputStream methods
			int numChannels = ais.getFormat().getChannels();
			int sampleSize = ais.getFormat().getSampleSizeInBits();
			boolean isBigEndian = ais.getFormat().isBigEndian();
			float sampleRate = ais.getFormat().getSampleRate();
			float frameRate = ais.getFormat().getFrameRate();
			int frameSize = ais.getFormat().getFrameSize();
			int frameLength = (int) ais.getFrameLength();

			// sampleRate = framerate = 44100.0 Hz (playback rate = sampling rate!)
			// 1 sec = 1000 millisecs
			// calculate delay in frames
			int frameDelay = (int) (timeDelay / 1000.0 * frameRate);

			// reset the AudioInputStream (mark goes to the start)
			ais.reset();

			// create a byte array of the right size
			// recall the lecture OHP slides ..
			a = new byte[(int) frameLength * frameSize];

			// fill the byte array with the data of the AudioInputStream
			ais.read(a);

			// Create an integer array, data, of the right size
			// only reason to do this is enabling type double mixing calculations
			// Each (channel) sample is made of 2 = sampleSize/8 bytes
			data = new int[a.length / 2];

			// fill the integer array by combining two 2 = sampleSize/8 bytes per sample of the
			// byte array a into one integer
			// Bytes HB and LB Big Endian make up one integer
			for (int i = 0; i < data.length; ++i) {
				/* First byte is HB (most significant digits) - coerce to 32-bit int */
				// HB =def sign_extend(a[2*i]) from 8 bit byte to 32 bit int
				int HB = (int) a[2 * i];
				/* Second byte is LB (least significant digits) - coerce to 32-bit int */
				// LB =def sign_extend(a[2*i+1]) from 8 bit byte to 32 bit int
				int LB = (int) a[2 * i + 1];
				// note that data[i] =def sign_extend(HB.LB)
				// | : Bool^32 x Bool^32 -----> Bool^32 where Bool = {0, 1}
				data[i] = HB << 8 | (LB & 0xff);
			}

			// split integer data array samples into two channels
			// if both channels are faded by the same factor
			// then there is no need to split the channels
			ch0 = new int[data.length / 2];
			ch1 = new int[data.length / 2];
			for (int i = 0; i < data.length / 2; i++) {
				ch0[i] = data[2 * i];
				ch1[i] = data[2 * i + 1];
			}

			// Adding a faded copy of the early signal to the later signal
			// THIS IS THE ECHO !!
			for (int i = frameDelay; i < ch0.length; ++i) {
				ch0[i] += (int) (ch0[i - frameDelay] * fading0);
				ch1[i] += (int) (ch1[i - frameDelay] * fading1);
			}

			// combine the two channels
			for (int i = 0; i < data.length; i += 2) {
				data[i] = ch0[i / 2];
				data[i + 1] = ch1[i / 2];
			}

			// get the maximum amplitute
			max = 0;
			for (int i = 0; i < data.length; ++i) {
				max = Math.max(max, Math.abs(data[i]));
			}

			// 16 digit 2s-complement range from -2^15 to +2^15-1 = 256*128-1
			// therefore we linearly scale data[i] values to lie within this range ..
			// .. so that each data[i] has a 16 digit "HB.LB binary representation"
			if (max > 256 * 128 - 1) {
				for (int i = 0; i < data.length; ++i) {
					data[i] = (int) (data[i] * (256.0 * 128.0 - 1) / max);
				}
			}

			// convert the integer array to a byte array
			for (int i = 0; i < data.length; ++i) {
				a[2 * i] = (byte) ((data[i] >> 8) & 0xff);
				a[2 * i + 1] = (byte) (data[i] & 0xff);
			}

		} catch (Exception e) {
			System.out.println("Something went wrong");
			e.printStackTrace();
		}

		// create a new AudioInputStream out of the the byteArray
		// and return it.
		return new AudioInputStream(new ByteArrayInputStream(a),
				ais.getFormat(), ais.getFrameLength());
	}

	/**** scaleToZero *****************************************************************/

	public static AudioInputStream scaleToZero(AudioInputStream ais) {


		byte[] a = null;
		int[] data;
		int[] ch0, ch1;

		int max;

		try {

			int frameSize = ais.getFormat().getFrameSize();
			int frameLength = (int) ais.getFrameLength();

			// reset the AudioInputStream (mark goes to the start) ??
			ais.reset();

			// create a byte array of the right size
			// recall the lecture OHP slides ..

			a = new byte[(int) frameLength * frameSize];

			// fill the byte array with the data of the AudioInputStream ??
			ais.read(a);

			// Create an integer array, data, of the right size
			// only reason to do this
			// is enabling type float/double mixing calculations ??
			data = new int[a.length / 2];

			// fill the integer array by combining two bytes of the
			// byte array a into
			// one integer - see lectures ??
			for (int i = 0; i < data.length; i++) {
				// same as echo
				int HB = (int) a[2 * i];
				int LB = (int) a[2 * i + 1];
				data[i] = HB << 8 | (LB & 0xff);
			}

			for (int i = 0; i < data.length / 2; i++) {

				double factor = 1 - ((3.0f / 4.0f) * (double) i / (double) (data.length / 2 - 1));
				data[2 * i] = (int) ((double) data[2 * i] * factor);
				data[2 * i + 1] = (int) ((double) data[2 * i + 1] * factor);
			}


			for (int i = 0; i < data.length; ++i) {
				a[2 * i] = (byte) ((data[i] >> 8) & 0xff);
				a[2 * i + 1] = (byte) (data[i] & 0xff);
			}

		} catch (Exception e) {
			System.out.println("Something went wrong");
			e.printStackTrace();
		}

		return new AudioInputStream(new ByteArrayInputStream(a), ais.getFormat(), ais.getFrameLength());

	}

	/**** addNote *****************************************************************/

	public static AudioInputStream addNote(AudioInputStream ais,
										   double frequency,
										   int noteLengthInMilliseconds) {
		byte[] a = null;
		int[] data;


		try {

			// number of frames for the note of noteLengthInMilliseconds
			float frameRate = ais.getFormat().getFrameRate();
			int frameSize = ais.getFormat().getFrameSize();
			int noteLengthInFrames = (int) ((float) noteLengthInMilliseconds / 1000.0f * frameRate);
			int noteLengthInBytes = noteLengthInFrames * frameSize;
			int noteLengthInInts = noteLengthInBytes / 2;

			a = new byte[noteLengthInBytes];
			data = new int[noteLengthInInts];
            // create the note as a data array of integer samples
            // each sample value data[i] is calculated using
            // the time t at which data[i] is played
			for (int i = 0; i < noteLengthInInts; i += 2) {
				// what is the time to play one frame?
				// BEFORE "frame" data[i]data[i+1] plays, how many frames are there?
				// hence compute t in terms of i
				double frameTime = 1 / frameRate;
				double doneFrames = i * 2 / frameSize;
				double t = frameTime * doneFrames;

				double amplitude = 5000;

				data[i] = (int) (amplitude * Math.sin(frequency * 2 * Math.PI * t));
				data[i + 1] = (int) (amplitude * Math.sin(frequency * 2 * Math.PI * t));

			}
            // copy the int data[i] array into byte a[i] array
			for (int i = 0; i < data.length; i++) {
				a[2 * i] = (byte) ((data[i] >> 8) & 0xff);
				a[2 * i + 1] = (byte) (data[i] & 0xff);
			}

		} catch (Exception e) {
			System.out.println("Something went wrong");
			e.printStackTrace();
		}


		return append(new AudioInputStream(new ByteArrayInputStream(a),
				ais.getFormat(), a.length / ais.getFormat().getFrameSize()), ais);

	}  // end addNote


	/**** append *****************************************************************/

	// THIS METHOD append IS SUPPLIED FOR YOU
	public static AudioInputStream append(AudioInputStream ais1, AudioInputStream ais2) {

		byte[] a, b, c = null;
		try {
			a = new byte[(int) ais1.getFrameLength() *
					ais1.getFormat().getFrameSize()];

			// fill the byte array with the data of the AudioInputStream
			ais1.read(a);
			b = new byte[(int) ais2.getFrameLength() *
					ais2.getFormat().getFrameSize()];

			// fill the byte array with the data of the AudioInputStream
			ais2.read(b);

			c = new byte[a.length + b.length];
			for (int i = 0; i < c.length; i++) {
				if (i < a.length) {
					c[i] = a[i];
				} else
					c[i] = b[i - a.length];
			}

		} catch (Exception e) {
			System.out.println("Something went wrong");
			e.printStackTrace();
		}


		return new AudioInputStream(new ByteArrayInputStream(c),
				ais1.getFormat(), c.length / ais1.getFormat().getFrameSize());
	} // end append

	/**** tune  *****************************************************************/

	public static AudioInputStream tune(AudioInputStream ais) {

		AudioInputStream temp = null;


		// create an empty AudioInputStream (of frame size 0)
		// call it temp (already declared above) 
		byte[] c = new byte[1];
		temp = new AudioInputStream(new ByteArrayInputStream(c), ais.getFormat(), 0);

		// specify variable names for both the frequencies in Hz and note lengths in seconds 
		// eg double C4, D4 etc for frequencies and s, l, ll, lll for lengths 
		// Hint: Each octave results in a doubling of frequency.

		double C4 = 261.63, E4 = 329.63, F4 = 349.23, G4 = 392.00, A4 = 440.00, B4 = 493.88, C5 = 523.25,
                D5 = 587.33, E5 = 2 * E4, Eb5 = 622.25, F5 = 2 * F4, G5 = 2 * G4, A5 = 2 * A4, B5 = 2 * B4,
                C6 = 2 * C5, D6 = 2 * D5, E6 = 2 * E5, F6 = 2 * F5, G6 = 2 * G5, A6 = 2 * A5, B6 = 2 * B5, C7 = 2 * C6;


		// and the lengths in milliseconds
		int s = 500, l = 2000, ll = 2500, lll = 2800;

		// also sprach zarathustra: 2001 A Space Odyssey 
		// specify the tune
		double[][] notes = {
				{C4, l}, {G4, l}, {C5, l},
				{E5, s}, {Eb5, lll},
				{C4, l}, {G4, l}, {C5, l},
				{Eb5, s}, {E5, lll},
				{A5, s}, {B5, s}, {C6, l},
				{A5, s}, {B5, s}, {C6, l},
				{D6, ll},
				{E6, s}, {F6, s}, {G6, l},
				{E6, s}, {F6, s}, {G6, l},
				{A6, l}, {B6, l}, {C7, lll}
		};

		// use addNote to build the tune as an AudioInputStream
		// starting from the empty AudioInputStream temp (above) and adding each note one by one using A LOOP 
		for (int i = notes.length - 1; i >= 0; i--) {
			temp = addNote(temp, notes[i][0], (int) notes[i][1]);
			if (i == 5) {
                temp = addNote(temp, 0, 500);
            }else if ( i == 10){
                temp = addNote(temp, 0, 500);
            } else {
				temp = addNote(temp, 0, 100);}
		}


		// append temp, ie the tune, to current ais 
		return append(temp, ais);
	}

	/**** altChannels *****************************************************************/

	public static AudioInputStream altChannels(AudioInputStream ais, double timeInterval) {


		int frameSize = ais.getFormat().getFrameSize(); // = 4
		float frameRate = ais.getFormat().getFrameRate();
		int inputLengthInFrames = (int) ais.getFrameLength();
		int frameInterval = (int) (frameRate * timeInterval); // number of frames played during timeInterval
		int inputLengthInBytes = frameSize * inputLengthInFrames;
		int numChannels = ais.getFormat().getChannels(); // = 2
        int N = frameInterval * frameSize / 2;

		// byte arrays for input channels and output channels
        int[] ich0, ich1, och0, och1, data = null;
        byte[] a, b = null;
		try{

            // create new byte arrays a for input and b for output of the right size
			a = new byte[(int) inputLengthInBytes];
			b = new byte[(int) inputLengthInBytes * 2];

            // fill the byte array a with the data of the AudioInputStream
			ais.read(a);

			data = new int[a.length / 2];
            // create new byte arrays for input and output channels of the right size
            // fill up ich0 and ich1 by splitting a
            for (int i = 0; i < data.length; i++){
				int HB = (int) a[2 * i];
				int LB = (int) a[2 * i + 1];
				data[i] = HB << 8 | (LB & 0xff);
			}


			ich0 = new int[data.length / numChannels];
			ich1 = new int[data.length / numChannels];

			och0 = new int[ich0.length * 2];
			och1 = new int[och0.length];

			int L = ich0.length;
			int outL = L * 2;
			int R = (outL % (2 * N)) / 2;

            for (int i = 0; i < data.length / 2; i++){
				ich0[i] = data[i * 2];
				ich1[i] = data[i * 2 + 1];
			}

			for (int i = 0; i < (double) outL / (double) (2 * N) - 1; i++) {
				for (int j = 0; j < N; j++) {
					och0[(i * 2 * N) + j] = ich0[(i * N) + j];
					och1[(i * 2 * N) + j] = 0x00;
					och0[(i * 2 * N) + j + N] = 0x00;
					och1[(i * 2 * N) + j + N] = ich1[(i * N) + j];
				}
			}

			int lastSeg = (int) Math.floor((double) outL / (double) (2 * N));
			for (int k = 0; k < R; k++){
				och0[(lastSeg * 2 * N) + k] = ich0[(lastSeg * N) + k];
				och1[(lastSeg * 2 * N) + k] = 0x00;

				och0[(lastSeg * 2 * N) + k + R] = 0x00;
				och1[(lastSeg * 2 * N) + k + R] = ich1[(lastSeg * N) + k];
			}

            // fill up b using och0 and och1
			data = new int[och0.length * 2];
			for (int i = 0; i < data.length; i += 2){
				data[i] = och0[i / 2];
				data[i + 1] = och1[i / 2];
			}

			// convert the integer array to a byte array
			for (int i = 0; i < data.length; i++) {
                    b[2 * i] = (byte) ((data[i] >> 8) & 0xff);
                    b[2 * i + 1] = (byte) (data[i] & 0xff);


            }

		}catch (Exception e){
            System.out.println("Something went wrong");
            e.printStackTrace();
        }

        // return b
        return new AudioInputStream(new ByteArrayInputStream(b), ais.getFormat(),
                b.length / ais.getFormat().getFrameSize());

    } // end altChannels

} // AudioManipulation
