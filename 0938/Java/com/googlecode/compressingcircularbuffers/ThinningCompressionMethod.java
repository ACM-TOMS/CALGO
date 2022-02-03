package com.googlecode.compressingcircularbuffers;
/** Compresses by simply eliminating intermediary values. */
public class ThinningCompressionMethod implements CompressionMethod {
  int count = 0;
  double compressedValue = Double.NaN;
  public void update(double value) {
    if (0 == count) compressedValue = value;
    count++;
  }
  public void reset() { compressedValue = Double.NaN; count = 0;}
  public double getCompressedValue() { return compressedValue;}
  public double mergeValues(double older, double newer) {return older;}
  public int getCount() { return count; }
}
