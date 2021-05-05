#ifdef __CINT__
#pragma link C++ class sf_g::waveform;
#pragma link C++ class sf_g::amplitude;
#pragma link C++ class sf_g::baseline;
#pragma link C++ class sf_g::cfd_time;
#pragma link C++ class sf_g::charge;
#pragma link C++ class sf_g::rise_time;
#pragma link C++ class sf_g::composite< sf_g::amplitude, sf_g::charge, sf_g::baseline >;
#pragma link C++ class sf_g::composite< sf_g::amplitude, sf_g::baseline >;
#pragma link C++ class sf_g::composite< sf_g::rise_time >;

#endif
