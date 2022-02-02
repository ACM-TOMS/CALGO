package com.googlecode.compressingcircularbuffers;

// tests the compressing circular buffer package.
public class CompressingCircularBufferTester  {

  private static void assertEq(double x1, double x2) {
    final double RELATIVE_DIFFERENCE = 0.01;
    if (Double.isNaN(x1) || Double.isNaN(x2))  {
      if (Double.isNaN(x1) != Double.isNaN(x2))
        throw new IllegalStateException("x1 ("+x1+") != x2 (" + x2 + ")");       
    }
    else if (Math.abs(x1 - x2) > RELATIVE_DIFFERENCE*(Math.abs(x1)+Math.abs(x2)))
      throw new IllegalStateException("x1 ("+x1+") != x2 (" + x2 + ")");
  }
  /**
   * Compares performance of CompressingCircularBuffers to
   * SimpleCompressingBuffers.
   * <p>
   *
   * This code was used to generate the performance comparison in the
   * article; it calculates the O(n) worst case update times of the
   * simple method, and validates that average case update times are
   * O(1) (and about the same) for both methods--used to check that the
   * CCB method does not add an unacceptably large fixed cost due to its
   * more complicated indexing.
   *
   * Note: you may need to add a java command line (VM) argument to
   * increase the heap such as 'java -Xmx1000000000'
   * if you get an error about not enough heap space when running this test.
   *
   */
  public static void timingTest() {
    final int COMPRESSIONS = 4;
    // number of multiples of the buffer size that the data sequence is
    final int ROUNDS = (int) Math.pow(2, COMPRESSIONS-1);
    final int[] BUFFER_SIZES = {1000001, 2000001, 4000001, 8000001, 16000001, };
    final int[] REPS =         {16,      8,       4,       2,        1};
    System.out.println("BufferSize, CircularAvg, CircularWorst, SimpleAvg, SimpleWorst"); 
    for (int i = 0; i < BUFFER_SIZES.length; i++) {
      CompressingCircularBuffer cb = new CompressingCircularBuffer(BUFFER_SIZES[i]);
      SimpleCompressingBuffer sb = new SimpleCompressingBuffer(BUFFER_SIZES[i]);
      for (int isSimple = 0; isSimple < 2; isSimple++) { 
        double t0Total = 0;
        double t1Total = 0;
        long t0 = System.currentTimeMillis();
        for (int j = 0; j < REPS[i]; j++) {
          int m = 1;  // helps determine when O(n) step will kick in for simple method
          if (isSimple == 1) {
            sb.reset();
          }
          else {
            cb.reset();
          }
          // the extra aggregate ("+ROUNDS") is needed to push out the
          // last compression sweep for the simple method.
          for (int k=0; k < ROUNDS*BUFFER_SIZES[i]+ROUNDS; k++) {
            if (k == m*BUFFER_SIZES[i]+m-1) { // O(n) step for simple method
              long t1 = System.currentTimeMillis();  
              if (isSimple == 1)
                sb.update(k);
              else
                cb.update(k);
              m *= 2;
              t1Total += System.currentTimeMillis() - t1;
            }
            else {
              if (isSimple == 1)
                sb.update(k);
              else
                cb.update(k);
            }
          }
        }
        t0Total = System.currentTimeMillis() - t0;
        if (isSimple == 0)
          System.out.print(BUFFER_SIZES[i] + ", " +
                           t0Total/(BUFFER_SIZES[i]*ROUNDS*REPS[i]) + ", " +
                           t1Total/(COMPRESSIONS*REPS[i]));
        else       
          System.out.println(", " + t0Total/(BUFFER_SIZES[i]*ROUNDS*REPS[i]) + ", " +
                             t1Total/(COMPRESSIONS*REPS[i]));
      }
    }
  }
  /**
   * Tests the CompressingCircularBuffer facility, and then runs 
   * the timingTest method.
   */
  public static void main(String[] args) {
    
   double[] values = {0, 1, 2};
   int nSamples = 0;
   // test with smallest non-trivial buffer (N=3).
	 CompressingCircularBuffer ccb = new CompressingCircularBuffer(3);
     for (int i = 0; i < 3; i++) {
       assertEq(ccb.getNCompressed(), i);
       assertEq(ccb.getNSamples(), nSamples);
	     ccb.update(i); nSamples++;
       assertEq(ccb.getCompressedValue(i),values[i]);
       assertEq(ccb.getCompressedIndex(i),i);
       assertEq(ccb.getFirstSampleIndex(i), i);
       assertEq(ccb.getLastSampleIndex(i), i);
     } 
     assertEq(ccb.getNCompressed(), 3);
     assertEq(ccb.getCompressionRatio(), 1);
     ccb.update(3); nSamples++;
     assertEq(ccb.getFirstSampleIndex(2),2);
     assertEq(ccb.getLastSampleIndex(2), 3);
     assertEq(ccb.getNSamples(), nSamples);
     assertEq(ccb.getNCompressedSamples(), nSamples);
     values[2] = (2 + 3)/2.;
     for (int i = 0; i < 3; i++) {
        assertEq(ccb.getCompressedValue(i), values[i]);
        assertEq(ccb.getCompressedIndex(i),i);
        assertEq(ccb.getNCompressed(), 3);
      }
     assertEq(ccb.getCompressionRatio(), 2);
     ccb.update(4); nSamples++;
     assertEq(ccb.getNSamples(), nSamples);
     assertEq(ccb.getNCompressedSamples(), nSamples-1);
     assertEq(ccb.getIncompleteCompressedValue(),4);
     ccb.update(5); nSamples++;
     assertEq(ccb.getFirstSampleIndex(2), 4);
     assertEq(ccb.getLastSampleIndex(2), 5);
     assertEq(ccb.getNSamples(), nSamples);
     assertEq(ccb.getNCompressedSamples(), nSamples);
     values[0] = (0 + 1)/2.;
     values[1] = (2 + 3)/2.;
     values[2] = (4 + 5)/2.;
     for (int i = 0; i < 3; i++) {
       assertEq(ccb.getCompressedValue(i),values[i]);
       assertEq(ccb.getCompressedIndex(2*i),i);
       assertEq(ccb.getCompressedIndex(2*i+1),i);
     }
     assertEq(ccb.getCompressionRatio(), 2);     

   // Something simple, kind of useful, and w lots of numbers in the
   // sequence: average, stdDev of a uniform distribution.  At end, each average
   // will have 65536 samples: should reliably get average, stdDev within 1%
	 CompressingCircularBuffer avg = new CompressingCircularBuffer(3);
	 CompressingCircularBuffer avg2 = new CompressingCircularBuffer(3);
	 for (int i = 0; i < 3*65536; i++) {
		double x = 2*Math.random();
		avg.update(x); 
		avg2.update(x*x); 
	 }
	
	 for (int i = 0; i < 3; i++) {
		 double avgi = avg.getCompressedValue(i);
		 double avgSquares = avg2.getCompressedValue(i);
		 assertEq(avgi, 1.0);
		 // Basic calculus: uniform dist 0 to 2 has variance of 1/3. Use well known
		 // stat formula: "variance = square of avg minus avg of squares" to est. variance.
		 assertEq(avgSquares-avgi*avgi, 1/3.);
   }

   // series of simple tests to exercise each built-in compression method
   ccb = new CompressingCircularBuffer(3, new ThinningCompressionMethod()); 
   for (int i = 0; i < 4*3; i++) {
     ccb.update(i);
   }
   
   for (int i = 0; i < 3; i++) {
     assertEq(ccb.getCompressedValue(i), i*4);
   }
   
   ccb = new CompressingCircularBuffer(3, new SummingCompressionMethod()); 
   for (int i = 0; i < 4*3; i++) {
     ccb.update(1);
   }
   for (int i = 0; i < 3; i++) 
     assertEq(ccb.getCompressedValue(i), 4);
   
   ccb = new CompressingCircularBuffer(3, new MinimizingCompressionMethod()); 
   for (int i = 0; i < 4*3; i++) {
     ccb.update(i);
   }
   for (int i = 0; i < 3; i++) 
     assertEq(ccb.getCompressedValue(i), i*4);
   ccb = new CompressingCircularBuffer(3, new MaximizingCompressionMethod()); 
   for (int i = 0; i < 4*3; i++) {
     ccb.update(i);
   }
   for (int i = 0; i < 3; i++) 
     assertEq(ccb.getCompressedValue(i), i*4+3);

   ccb = new CompressingCircularBuffer(65535, new RandomSamplingCompressionMethod()); 
   for (int i = 0; i < 65535*8; i++) {
     ccb.update(i%8);
   }
   double sum = 0;
   for (int i = 0; i < 65535; i++)  
     sum += ccb.getCompressedValue(i);

   assertEq(sum/65535., 3.5);


   // test to exercise min and max compression ratio features
   ccb = new CompressingCircularBuffer(5, new SummingCompressionMethod(),
                                       3, 6);
   for (int i = 0; i < 15; i++)
     ccb.update(1);
   
   for (int i = 0; i < 5; i++) {
     assertEq(ccb.getCompressedValue(i), 3);
     assertEq(ccb.getCompressedIndex(3*i),i);
     assertEq(ccb.getCompressedIndex(3*i+1),i);
     assertEq(ccb.getCompressedIndex(3*i+2),i);
   }
   for (int i = 0; i < 15; i++)
     ccb.update(1);

   for (int i = 0; i < 5; i++) {
     assertEq(ccb.getCompressedValue(i), 6);
     assertEq(ccb.getCompressedIndex(6*i),i);
     assertEq(ccb.getCompressedIndex(6*i+1),i);
     assertEq(ccb.getCompressedIndex(6*i+5),i);
   }
   
   for (int i = 0; i < 6; i++)
     ccb.update(1);

   assertEq(ccb.getCompressedIndex(0),CompressingCircularBuffer.NOSUCHINDEX);
   assertEq(ccb.getCompressedIndex(5),CompressingCircularBuffer.NOSUCHINDEX);
   for (int i = 0; i < 5; i++) {
     assertEq(ccb.getCompressedValue(i), 6);
     assertEq(ccb.getCompressedIndex(6*(i+1)),i);
     assertEq(ccb.getCompressedIndex(6*(i+1)+1),i);
     assertEq(ccb.getCompressedIndex(6*(i+1)+5),i);
   }
   
   for (int i = 0; i < 4*6; i++)
     ccb.update(1);

   assertEq(ccb.getCompressedIndex(0),CompressingCircularBuffer.NOSUCHINDEX);
   assertEq(ccb.getCompressedIndex(4*6+5),CompressingCircularBuffer.NOSUCHINDEX);
   for (int i = 0; i < 5; i++) {
     assertEq(ccb.getCompressedValue(i), 6);
     assertEq(ccb.getCompressedIndex(6*(i+5)),i);
     assertEq(ccb.getCompressedIndex(6*(i+5)+1),i);
     assertEq(ccb.getCompressedIndex(6*(i+5)+5),i);
   }
   // Test that all variants of constructors work consistently w
   // default parameters. Also tests the special case of a 1-element
   // sized buffer.

   CompressingCircularBuffer[] ccbArray = {
     new CompressingCircularBuffer(),
     new CompressingCircularBuffer(1),
     new CompressingCircularBuffer(1, new AveragingCompressionMethod()),
     new CompressingCircularBuffer(1, new AveragingCompressionMethod(), 1),
     new CompressingCircularBuffer(1, new AveragingCompressionMethod(), 1, CompressingCircularBuffer.DEFAULT_MAX_COMPRESSION_RATIO)
   };

   sum = 0;
   int count = 0;
   for (int j = 0; j < ccbArray.length; j++)
     ccbArray[j].update(0);
   count++;
   for (int k = 0; k < 8; k++) {
     for (int i = 0; i < Math.pow(2,k); i++) {
       for (int j = 0; j < ccbArray.length; j++)
         ccbArray[j].update(i);
       sum += i;
       count++;
     }
     for (int j = 0; j < ccbArray.length; j++) {
       assertEq(ccbArray[j].getCompressedValue(0),sum/count);
       assertEq(ccbArray[j].getCompressionRatio(),  Math.pow(2,k+1));
       assertEq(ccbArray[j].getFirstSampleIndex(0),0);
       assertEq(ccbArray[j].getLastSampleIndex(0),Math.pow(2,k+1)-1);
       assertEq(ccbArray[j].getIncompleteCompressedValue(),Double.NaN);       
       assertEq(ccbArray[j].getNSamples(),Math.pow(2,k+1));       
       assertEq(ccbArray[j].getNCompressedSamples(),Math.pow(2,k+1));       
       assertEq(ccbArray[j].getNCompressed(),1);
     }
   }

   // test with minCompressionRatio > 1, and what happens when max
   // compression ratio is exceeded

   ccb = new CompressingCircularBuffer(5, new SummingCompressionMethod(), 5, 5);
   for (int i = 0; i < 10; i++)
     for (int k = 0; k < 5; k++)
       ccb.update(i-2+k);         // average value is i

   for (int i = 5; i < 10; i++) 
     assertEq(ccb.getCompressedValue(i-5), i*5);

// generate relative performance data quoted in article.
   timingTest();
  }
}
