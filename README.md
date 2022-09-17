## British National Grid (BNG) to Latitude and Longitude
This is a straight up port - comments and all - of the [bng_to_latlon Python library](https://github.com/fmalina/blocl-bnglatlon) by [Hannah Fry](https://hannahfry.co.uk/) and F. Malina.

>Converts British National Grid eastings and northings (OSBG36) to latitude and longitude (WGS84) ~~and vice versa, as used by https://blocl.uk~~.

## Usage
```golang
package main

import (
	"fmt"

	"github.com/jkfsom/bnglatlng"
)

func main(){
	easting, northing := 382839.0, 150851.36084

	lat, lng := bnglatlng.OSGB36toWGS84(easting, northing)
	fmt.Printf("easting: %f, northing: %f -> latitude: %f, longitude: %f\n", easting, northing, lat, lng)
}
```