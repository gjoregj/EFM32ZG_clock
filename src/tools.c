/*
 * tools.c
 *
 *  Created on: 28.2.2014
 *      Author: GJOREGJ
 */

#include <time.h>
#include <math.h>

#define ONE_DAY 	86400
#define TROP_YEAR 	31556925
#define ANOM_YEAR 	31558433
#define INCLINATION 0.409105176667471		/* Earths axial tilt at the epoch */
#define PERIHELION 	31316400    			/* perihelion of 1999, 03 jan 13:00 UTC */
#define SOLSTICE 	836160        			/* winter solstice of 1999, 22 Dec 07:44 UTC */
#define TWO_PI 		6.283185307179586
#define TROP_CYCLE 	5022440.6025
#define ANOM_CYCLE 	5022680.6082
#define DELTA_V 	0.03342044   			/* 2x orbital eccentricity */

#define LAG			38520

static volatile long __latitude = 42.0;
static volatile long __longitude = 21.4333;

int equation_of_time(const time_t * timer)
{
  int s, p;
  double pf, sf, dV;

  /* compute orbital position relative to perihelion */
  p = *timer % ANOM_YEAR;
  p += PERIHELION;
  pf = p;
  pf /= ANOM_CYCLE;
  pf = sin(pf);

  /* Derive a velocity correction factor from the perihelion angle */
  dV = pf * DELTA_V;

  /* compute approximate position relative to solstice */
  s = *timer % TROP_YEAR;
  s += SOLSTICE;
  s *= 2;
  sf = s;
  sf /= TROP_CYCLE;

  /* modulate to derive actual position */
  sf += dV;
  sf = sin(sf);

  /* compute contributions */
  sf *= 592.2;
  pf *= 459.6;
  s = pf + sf;
  return -s;

}

time_t solar_noon(const time_t * timer)
{
  time_t t;
  long n;

  /* determine time of solar noon at the prime meridian */
  t = *timer % ONE_DAY;
  t = *timer - t;
  t += 43200L;
  t -= equation_of_time(timer);

  /* rotate to observers longitude */
  n = __longitude / 15L;
  t -= n;

  return t;

}

double solar_declination(const time_t * timer)
{

  unsigned int fT, oV;
  double dV, dT;

  /* Determine orbital angle relative to perihelion of January 1999 */
  oV = *timer % ANOM_YEAR;
  oV += PERIHELION;
  dV = oV;
  dV /= ANOM_CYCLE;

  /* Derive velocity correction factor from the perihelion angle */
  dV = sin(dV);
  dV *= DELTA_V;

  /* Determine orbital angle relative to solstice of December 1999 */
  fT = *timer % TROP_YEAR;
  fT += SOLSTICE + LAG;
  dT = fT;
  dT /= TROP_CYCLE;
  dT += dV;

  /* Finally having the solstice angle, we can compute the declination */
  dT = cos(dT) * INCLINATION;

  return -dT;
}

long daylight_seconds(const time_t * timer)
{
  double l, d;
  long n;

  /* convert latitude to radians */
  l = __latitude / 206264.806;

  d = -solar_declination(timer);

  /* partial 'Sunrise Equation' */
  d = tan(l) * tan(d);

  /* magnitude of d may exceed 1.0 at near solstices */
  if (d > 1.0)
    d = 1.0;

  if (d < -1.0)
    d = -1.0;

  /* derive hour angle */
  d = acos(d);

  /* but for atmospheric refraction, this would be d /= M_PI */
  d /= 3.112505;

  n = ONE_DAY * d;

  return n;
}

time_t sun_rise(const time_t * timer)
{
  long n;
  time_t t;

  /* sunrise is 1/2 'day' before solar noon */
  t = solar_noon(timer);
  n = daylight_seconds(timer) / 2L;
  t -= n;

  return t;
}

time_t sun_set(const time_t * timer)
{
  long n;
  time_t t;

  /* sunset is 1/2 'day' after solar noon */
  t = solar_noon(timer);
  n = daylight_seconds(timer) / 2L;
  t += n;

  return t;

}

void set_position(long lat, long lon)
{
  __latitude = lat;
  __longitude = lon;
}
