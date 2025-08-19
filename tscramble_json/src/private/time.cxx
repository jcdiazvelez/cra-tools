#include "astro/time.h"
#include "units.h"
#include "star/pal.h"

#include <string>
#include <cstdio>


// Valid range check //

bool astro::Time::IsValid() {
  if (fractionMJD_ < 0. || fractionMJD_ >= 1.) {
    printf("fractionMJD_ outside range [0.,1.)\n");
    return false;
  }
  return true;
}


// Set Time Methods //

void astro::Time::SetTime(double mjd) {
  // n.b. use floor() for proper handling if mjd is negative
  wholeMJD_ = int(floor(mjd)); // example: int(floor(-2.5)) = -3 
  fractionMJD_ = mjd - wholeMJD_;
  IsValid();
}

void astro::Time::SetTime(int wholeMJD, double fractionMJD) {
  wholeMJD_ = wholeMJD;
  fractionMJD_ = fractionMJD;
  IsValid();
}
      
void astro::Time::SetTime(int MJday, int MJsec, double MJns) {
  wholeMJD_ = MJday;
  fractionMJD_ = (MJsec + 1.e-9*MJns)/86400.;
  IsValid();
}

void astro::Time::SetTime(int year, int month, int day, 
			  int hour, int min, double sec) {
  int status;
  double wholeMJD, fractionMJD;
  
  palCldj(year, month, day, &wholeMJD, &status);
  if (status) 
    { printf("Invalid year/month/day in SetTime\n"); }
  
  palDtf2d(hour, min, sec, &fractionMJD, &status);
  if (status) 
    { printf("Invalid hour,min,sec in SetTime\n"); }
  
  SetTime(int(wholeMJD), fractionMJD);
}

void astro::Time::SetTime_JulianDay(double julianDay) { 
  SetTime(julianDay - astro::mjd_offset); 
}

void astro::Time::SetTime_JEpoch(double jEpoch) { 
  SetTime( palEpj2d(jEpoch) ); 
}


// Get Time Methods //

void astro::Time::GetCalendarTime(int *year, int *month, int *day, 
				  int *hour, int *minute, 
				  double *second) const {
  int status;
  double fractionDay;
  palDjcl( GetMJD(), year, month, day, &fractionDay, &status);
  *hour =   int( fractionMJD_*24. );
  *minute = int( fmod( fractionMJD_*24.*60. , 60.) );
  *second = fmod( fractionMJD_*86400. , 60.);
}

double astro::Time::GetTime_JulianDay() const { 
  return astro::mjd_offset + GetMJD(); 
}

double astro::Time::GetJEpoch() const { 
  return palEpj( GetMJD() ); 
}




// String Format Calendar Time

std::string astro::Time::GetCalendarString() {
  int year, month, day, hour, minute;
  double second;
  GetCalendarTime(&year, &month, &day, &hour, &minute, &second);
  char calFormat[100];
  sprintf(calFormat,"%4d/%02d/%02d %2d:%02d:%08.5f",
	  year, month, day, hour, minute, second);
  return std::string(calFormat);
}
