package bnglatlng

import (
	"fmt"
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestOSGB36toWGS84(t *testing.T) {
	cases := []struct {
		identifier  string
		easting     float64
		northing    float64
		expectedLat float64
		expectedLng float64
	}{
		{
			identifier:  "test case A",
			easting:     382839.0,
			northing:    150851.36084,
			expectedLat: 51.256564,
			expectedLng: -2.247300,
		},
	}

	for _, item := range cases {
		lat, lng := OSGB36toWGS84(item.easting, item.northing)
		assert.Equal(t, item.expectedLat, lat, fmt.Sprintf("test case ID '%s' - lat.", item.identifier))
		assert.Equal(t, item.expectedLng, lng, fmt.Sprintf("test case ID '%s' - lng.", item.identifier))
	}
}
