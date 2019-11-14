/**
 * @file maptiles.cpp
 * Code for the maptiles function.
 */

#include <iostream>
#include <map>
#include "maptiles.h"
//#include "cs225/RGB_HSL.h"

using namespace std;
using namespace cs225;

Point<3> convertToXYZ(LUVAPixel pixel) {
    return Point<3>( pixel.l, pixel.u, pixel.v );
}

MosaicCanvas* mapTiles(SourceImage const& theSource,
                       vector<TileImage>& theTiles)
{
    /**
     * @todo Implement this function!
     */
    vector<Point<3>> kdtree_grid;
    map<Point<3>, int> color_map;
    for (size_t i = 0; i < theTiles.size(); i++) {
      LUVAPixel color_tiles = theTiles[i].getAverageColor();
      Point<3> colorTiles = convertToXYZ(color_tiles);
      pair<Point<3>, TileImage> element(colorTiles,theTiles[i]);
      kdtree_grid.push_back(colorTiles);
      color_map[colorTiles] = i;
    }
    KDTree<3> targetTree(kdtree_grid);
    int row = theSource.getRows();
    int column = theSource.getColumns();
    MosaicCanvas* outputPtr = new MosaicCanvas(row, column);
    for (int x = 0; x < row; x++) {
      for (int y = 0; y < column; y++) {
        LUVAPixel pixel = theSource.getRegionColor(x,y);
        Point<3> point_pixel = convertToXYZ(pixel);
        Point<3> closest = targetTree.findNearestNeighbor(point_pixel);
        int targetTileImageNUM = color_map[closest];
        TileImage *targetTileImage = &(theTiles[targetTileImageNUM]);
        outputPtr->setTile(x,y,targetTileImage);
      }
    }
    return outputPtr;
}
