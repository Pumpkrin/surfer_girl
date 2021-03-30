surfer_girl is a tool suit dedicated to extraction of binary data produced by the wavecatcher system (https://owncloud.lal.in2p3.fr/index.php/s/56e4a2c53a991cb08f73d03f1ce58ba2?path=%2F).
It is composed, for now, of two different parts:
 - surfer_girl, which is used to transform the raw data into a ROOT(cern) compatible format, i.e.: histograms (TH1D).
   This executable takes a total of 4 arguments: you can specify which file to read with : -in input_file.bin, and where to output with -out output_file.root.
 - tree_flip, which should be used to extract information from the waveform in histogram form. 
   It will produce exactly what you ask of it, meaning that it will only extract the information you specify for each channel.
   It should be launched the following way: 
     ./tree_flip -in input_file1.root:input_file2.root:... -out output_file.root -mod {n:ms...}...
   Here several input files can be read and concatenated in the same output.
   The -mod argument is used to specify the information you want to extract from the files:
      it should contain the channel number from which you want to extract the data ( i.e. n above )
      and the list of modules you want to apply to it, which are for now: a(amplitude), b(baseline), c(charge), r(rise time), t(cfd time)
  
Its use requires c++14 (and a compliant ROOT installation). CMake should be used to build the executables.
