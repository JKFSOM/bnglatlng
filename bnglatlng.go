/*
This is a straight up port - comments and all - of the bng_to_latlon Python
library by Hannah Fry and F. Malina.
---
https://github.com/fmalina/blocl-bnglatlon
*/
package bnglatlng

import "math"

// convert seconds to radians
func sec_to_rad(x float64) float64 {
	return x * math.Pi / (180 * 3600.)
}

func OSGB36toWGS84(E, N float64) (float64, float64) {
	// The Airy 1830 semi-major and semi-minor axes used for OSGB36 (m)
	a, b := 6377563.396, 6356256.909
	F0 := 0.9996012717 // scale factor on the central meridian

	// Latitude and longtitude of true origin (radians)
	lat0 := 49 * math.Pi / 180
	lon0 := -2 * math.Pi / 180 // longtitude of central meridian

	// Northing & easting of true origin (m)
	N0, E0 := -100000.0, 400000.0
	e2 := 1 - (b*b)/(a*a) // eccentricity squared
	n := (a - b) / (a + b)

	// Initialise the iterative variables
	lat, M := lat0, 0.0

	for N-N0-M >= 0.00001 { //Accurate to 0.01mm
		lat = (N-N0-M)/(a*F0) + lat
		M1 := (1 + n + (5./4)*math.Pow(n, 2) + (5./4)*math.Pow(n, 3)) * (lat - lat0)
		M2 := (3*n + 3*math.Pow(n, 2) + (21./8)*math.Pow(n, 3)) * math.Sin(lat-lat0) * math.Cos(lat+lat0)
		M3 := ((15./8)*math.Pow(n, 2) + (15./8)*math.Pow(n, 3)) * math.Sin(2*(lat-lat0)) * math.Cos(2*(lat+lat0))
		M4 := (35. / 24) * math.Pow(n, 3) * math.Sin(3*(lat-lat0)) * math.Cos(3*(lat+lat0))
		// meridional arc
		M = b * F0 * (M1 - M2 + M3 - M4)
	}

	// transverse radius of curvature
	nu := a * F0 / math.Sqrt(1-e2*math.Pow(math.Sin(lat), 2))

	//meridional radius of curvature
	rho := a * F0 * (1 - e2) * math.Pow((1-e2*math.Pow(math.Sin(lat), 2)), -1.5)
	eta2 := nu/rho - 1

	sec_lat := 1. / math.Cos(lat)
	VII := math.Tan(lat) / (2 * rho * nu)
	VIII := math.Tan(lat) / (24 * rho * math.Pow(nu, 3)) * (5 + 3*math.Pow(math.Tan(lat), 2) + eta2 - 9*math.Pow(math.Tan(lat), 2*eta2))
	IX := math.Tan(lat) / (720 * rho * math.Pow(nu, 5)) * (61 + 90*math.Pow(math.Tan(lat), 2) + 45*math.Pow(math.Tan(lat), 4))
	X := sec_lat / nu
	XI := sec_lat / (6 * math.Pow(nu, 3)) * (nu/rho + 2*math.Pow(math.Tan(lat), 2))
	XII := sec_lat / (120 * math.Pow(nu, 5)) * (5 + 28*math.Pow(math.Tan(lat), 2) + 24*math.Pow(math.Tan(lat), 4))
	XIIA := sec_lat / (5040 * math.Pow(nu, 7)) * (61 + 662*math.Pow(math.Tan(lat), 2) + 1320*math.Pow(math.Tan(lat), 4) + 720*math.Pow(math.Tan(lat), 6))
	dE := E - E0

	// These are on the wrong ellipsoid currently: Airy 1830 (denoted by _1)
	lat_1 := lat - VII*math.Pow(dE, 2) + VIII*math.Pow(dE, 4.) - IX*math.Pow(dE, 6.)
	lon_1 := lon0 + X*dE - XI*math.Pow(dE, 3) + XII*math.Pow(dE, 5.) - XIIA*math.Pow(dE, 7.)

	// Want to convert to the GRS80 ellipsoid.
	// First convert to cartesian from spherical polar coordinates
	H := 0.0 // Third spherical coord.
	x_1 := (nu/F0 + H) * math.Cos(lat_1) * math.Cos(lon_1)
	y_1 := (nu/F0 + H) * math.Cos(lat_1) * math.Sin(lon_1)
	z_1 := ((1-e2)*nu/F0 + H) * math.Sin(lat_1)

	// Perform Helmut transform (to go between Airy 1830 (_1) and GRS80 (_2))
	s := -20.4894 * math.Pow(10, -6) // The scale factor -1
	// The translations along x, y, z axes respectively
	tx, ty, tz := 446.448, -125.157, +542.060
	// The rotations along x, y, z respectively (in seconds)
	rxs, rys, rzs := 0.1502, 0.2470, 0.8421

	rx := sec_to_rad(rxs)
	ry := sec_to_rad(rys)
	rz := sec_to_rad(rzs)

	x_2 := tx + (1+s)*x_1 + (-rz)*y_1 + (ry)*z_1
	y_2 := ty + (rz)*x_1 + (1+s)*y_1 + (-rx)*z_1
	z_2 := tz + (-ry)*x_1 + (rx)*y_1 + (1+s)*z_1

	// Back to spherical polar coordinates from cartesian
	// Need some of the characteristics of the new ellipsoid

	// The GSR80 semi-major and semi-minor axes used for WGS84(m)
	a_2, b_2 := 6378137.000, 6356752.3141
	e2_2 := 1 - (b_2*b_2)/(a_2*a_2) // The eccentricity of the GRS80 ellipsoid
	p := math.Sqrt(math.Pow(x_2, 2) + math.Pow(y_2, 2))

	// Lat is obtained by an iterative proceedure:
	lat = math.Atan2(z_2, (p * (1 - e2_2))) // Initial value
	latold := 2 * math.Pi

	var nu_2 float64
	for math.Abs(lat-latold) > math.Pow(10, -16) {
		lat, latold = latold, lat
		nu_2 = a_2 / math.Sqrt(1-e2_2*math.Pow(math.Sin(latold), 2))
		lat = math.Atan2(z_2+e2_2*nu_2*math.Sin(latold), p)
	}

	lon := math.Atan2(y_2, x_2)

	// Convert to degrees
	lat = lat * 180 / math.Pi
	lon = lon * 180 / math.Pi

	return round(lat, 6), round(lon, 6)
}

func round(x float64, y int) float64 {
	places := math.Pow(10, float64(y))
	return math.Round(x*places) / places
}
