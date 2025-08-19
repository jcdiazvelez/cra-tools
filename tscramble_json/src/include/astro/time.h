#ifndef ASTRO_TIME_H_INCLUDED
#define ASTRO_TIME_H_INCLUDED

#include <string>
#include <math.h>



namespace astro
{

  static const double mjd_offset      = 2400000.5;
  static const double j2000_JulianDay = 2451545.0;


  class Time {

   private:
    int wholeMJD_;
    double fractionMJD_;
    
   public:

    // Constructor

    Time() 
      { SetTime(0.); }

    // Overload constructor with common ways for setting time

    Time(double mjd) 
      { SetTime(mjd); }

    Time(int wholeMJD, double fractionMJD) 
      { SetTime(wholeMJD, fractionMJD); }

    Time(int MJday, int MJsec, double MJns)
      { SetTime(MJday, MJsec, MJns); }

    Time(int year, int month, int day, int hour, int min, double sec) 
      { SetTime(year, month, day, hour, min, sec); }


    // Valid range check //

    bool IsValid();


    // Set Time Methods //

    void SetTime(double mjd);
    void SetTime(int wholeMJD, double fractionMJD);
    void SetTime(int MJday, int MJsec, double MJns);
    void SetTime(int year, int month, int day, int hour, int min, double sec);
    void SetTime_JulianDay(double julianDay);
    void SetTime_JEpoch(double jEpoch);


    // Get Time Methods //

    double GetMJD() const { return wholeMJD_ + fractionMJD_; }

    int GetWholeMJD() const { return wholeMJD_; }
    double GetFractionMJD() const { return fractionMJD_; }

    int GetMJday_Part() const { return wholeMJD_; }
    int GetMJsec_Part() const { return int(fractionMJD_ * 86400.); }
    double GetMJns_Part() const 
      { return 1.e9 * fmod(fractionMJD_*86400., 1.); }

    void GetCalendarTime(int *year, int *month, int *day, 
			 int *hour, int *minute, double *second) const;

    double GetTime_JulianDay() const;

    double GetJEpoch() const;


    // String Format Calendar Time

    std::string GetCalendarString();
  };

}

#endif // ASTRO_TIME_H_INCLUDED
