#include "Document.h"

#include <math.h>		/* floor */
#include <iostream>		/* cout */
#include <fstream>		/* fstream */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <map>          /* map */

#include "CDT/CDT.h"	/* triangulation */

bool Document::is_init = false;

int zonestat::cols = 0;
int zonestat::rows = 0;
double zonestat::width = 0;
double zonestat::height = 0;
double zonestat::max_entropy = -1;
double zonestat::most_common_relative_entropy = -1;
ubyte zonestat::most_common_entropy_byte = 0;

double* CartesianToBarycentric(const int x, const int y, const Point R1, const Point R2, const Point R3)
{
	double* coords = new double[3]{ 0.0 };

	double detT = (R2.y - R3.y) * (R1.x - R3.x) + (R3.x - R2.x) * (R1.y - R3.y);

	coords[0] = ((R2.y - R3.y) * (x - R3.x) + (R3.x - R2.x) * (y - R3.y)) / detT;
	coords[1] = ((R3.y - R1.y) * (x - R3.x) + (R1.x - R3.x) * (y - R3.y)) / detT;
	coords[2] = 1.0 - coords[0] - coords[1];

	return coords;
}

void Document::ILDataToRGB()
{
	initPixels();
	ilBindImage(image_name);
	if (ilGetInteger(IL_IMAGE_FORMAT) == IL_LUMINANCE)	//Document is greyscale
	{
		initGScPixels();
		ILubyte* data = new ILubyte[width * height * 1];
		ilCopyPixels(0, 0, 0, width, height, 1, IL_LUMINANCE, IL_UNSIGNED_BYTE, data);
		for (int i = 0; i < width; i++)
		{
			for (int j = 0; j < height; j++)
			{
				gscpixels[i][j] = data[j * width + i];
			}
		}
	}

	ILubyte* data = new ILubyte[width * height * 4];
	ilCopyPixels(0, 0, 0, width, height, 1, IL_RGBA, IL_UNSIGNED_BYTE, data);
	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			ubyte rgba[4] = { 0 };
			for (int k = 0; k < 4; k++)
				rgba[k] = data[(j * width + i) * 4 + k];
			rgbpixels[i][j] = RGBA(rgba[0], rgba[1], rgba[2], rgba[3]);
		}
	}
	delete[] data;
}

ILubyte* Document::RGBToILData(bool whiten) const
{
	ilBindImage(image_name);
	ILubyte* data = new ILubyte[width * height * 4];

	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			RGBA current_pixel = rgbpixels[i][j];

			if (whiten && !gscpixels.empty() && gscpixels[i][j] == GScPixel::WHITE)
			{
				current_pixel = RGBA::WHITE();
				current_pixel.alpha = rgbpixels[i][j].alpha;
			}

			int byte = 0;
			data[(j * width + i) * 4 + byte++] = current_pixel.red;
			data[(j * width + i) * 4 + byte++] = current_pixel.green;
			data[(j * width + i) * 4 + byte++] = current_pixel.blue;
			data[(j * width + i) * 4 + byte] = current_pixel.alpha;
		}
	}

	return data;
}

ILubyte* Document::GScToILData() const
{
	ILubyte* data = new ILubyte[width * height * 1];

	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			data[j * width + i] = gscpixels[i][j];
		}
	}

	return data;
}

void Document::initPixels()
{
	rgbpixels.clear();
	rgbpixels.resize(width);
	for (int i = 0; i < width; i++)
		rgbpixels[i].resize(height);
}

void Document::initGScPixels()
{
	gscpixels.clear();
	gscpixels.resize(width);
	for (int i = 0; i < width; i++)
		gscpixels[i].resize(height);
}

void Document::initFrequency(frequencies& frequency)
{
	constexpr int size = GScPixel::WHITE + 1;
	frequency.resize(size);
}

Document::Document(int width, int height)
{
	this->width = width;
	this->height = height;
	image_name = 0;
}

Document::Document(std::string file)
{
	height = width = 0;
	image_name = 0;
	image_path = file;
	Read();
}

Document::Document(ILuint image_name, std::string file)
{
	image_name = image_name;
	ilBindImage(image_name);
	width = ilGetInteger(IL_IMAGE_WIDTH);
	height = ilGetInteger(IL_IMAGE_HEIGHT);
	image_path = file;
	ILDataToRGB();
}

Document::Document(const Document& that)
{
	image_name = that.image_name;
	width = that.width;
	height = that.height;
	rgbpixels = that.rgbpixels;
	gscpixels = that.gscpixels;
	is_background_removed = that.is_background_removed;
}

void Document::init()
{
	if (ilGetInteger(IL_VERSION_NUM) < IL_VERSION)
	{
		printf("DevIL version is different...exiting!\n");
		exit(1);
	}

	ilInit();
	ilEnable(IL_FILE_OVERWRITE);
	ilEnable(IL_ORIGIN_SET);
	ilOriginFunc(IL_ORIGIN_UPPER_LEFT);

	is_init = true;
}

void Document::ValidateDimensions(int& minWidth, int& minHeight, int& maxWidth, int& maxHeight) const
{
	if (minWidth < 0) minWidth = 0;
	if (minWidth > width - 1) minWidth = width - 1;

	if (minHeight < 0) minHeight = 0;
	if (minHeight > height - 1) minHeight = height - 1;

	if (maxWidth == 0 || maxWidth <= minWidth || maxWidth > width) maxWidth = width;
	if (maxHeight == 0 || maxHeight <= minHeight || maxHeight > height) maxHeight = height;
}

Document& Document::operator=(const Document& that)
{
	image_name = that.image_name;
	width = that.width;
	height = that.height;
	rgbpixels = that.rgbpixels;
	is_background_removed = that.is_background_removed;
	gscpixels = that.gscpixels;
	return *this;
}

std::vector<std::vector<RGBA>>& Document::GetRGBPixels()
{
	return rgbpixels;
}

std::vector<std::vector<ubyte>>& Document::GetGScPixels()
{
	return gscpixels;
}

void Document::RGBtoGSc()
{
	initGScPixels();
	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			gscpixels[i][j] = rgbpixels[i][j].ToGrey();
		}
	}
}

Document Document::Crop(int minWidth, int minHeight, int maxWidth, int maxHeight) const
{
	ValidateDimensions(minWidth, minHeight, maxWidth, maxHeight);

	int cropWidth = maxWidth - minWidth;
	int cropHeight = maxHeight - minHeight;

	ilBindImage(image_name);
	int type = ilGetInteger(IL_IMAGE_TYPE);
	int format = ilGetInteger(IL_IMAGE_FORMAT);

	ILuint new_name;
	ilGenImages(1, &new_name);
	ilBindImage(new_name);

	bool success = false;
	success = ilTexImage(cropWidth, cropHeight, 0, 4, format, type, NULL);
	ilRegisterOrigin(IL_ORIGIN_UPPER_LEFT);

	Document Crop = Document(new_name, image_path);
	std::vector<std::vector<RGBA>>& rgbpx = Crop.GetRGBPixels();
	std::vector<std::vector<ubyte>>& greypx = Crop.GetGScPixels();
	for (int i = 0; i < cropWidth; i++)
	{
		for (int j = 0; j < cropHeight; j++)
		{
			rgbpx[i][j] = rgbpixels[minWidth + i][minHeight + j];
		}
	}

	if (!gscpixels.empty())
	{
		for (int i = 0; i < cropWidth; i++)
		{
			for (int j = 0; j < cropHeight; j++)
				greypx[i][j] = gscpixels[minWidth + i][minHeight + j];
		}
	}

	return Crop;
}

frequencies Document::GetGScFrequency(int minWidth, int minHeight, int maxWidth, int maxHeight)
{
	if (gscpixels.empty()) RGBtoGSc();
	frequencies frequency;
	initFrequency(frequency);
	ValidateDimensions(minWidth, minHeight, maxWidth, maxHeight);
	for (int i = minWidth; i < maxWidth; i++)
	{
		for (int j = minHeight; j < maxHeight; j++)
		{
			unsigned char idx = gscpixels[i][j];
			frequency[idx]++;
		}
	}
	return frequency;
}

int Document::InitZones(int amount)
{
	if (width < MIN_IMAGE_WIDTH || height < MIN_IMAGE_HEIGHT)
		throw image_too_small(("Image resolution must be at least " + std::to_string(MIN_IMAGE_WIDTH) + "x" + std::to_string(MIN_IMAGE_HEIGHT) + "!").c_str());

	//Calculating zone borders:
	const unsigned long long int image_area = (long long)(width) * (long long)(height);
	double zone_area = (double)image_area / amount;

	double zone_side_length = sqrt(zone_area);
	if (zone_side_length < MIN_ZONE_SIDE_LENGTH || amount == 0)
		zone_side_length = MIN_ZONE_SIDE_LENGTH;
	if (MAX_ZONE_SIDE_LENGTH < zone_side_length)
		zone_side_length = MAX_ZONE_SIDE_LENGTH;

	const int cols = width / zone_side_length;
	const int rows = height / zone_side_length;

	amount = cols * rows;

	const double zone_width = double(width) / cols;
	const double zone_height = double(height) / rows;

	zonestat::cols = cols;
	zonestat::rows = rows;
	zonestat::width = zone_width;
	zonestat::height = zone_height;

	zones.resize(cols);
	for (size_t i = 0; i < zones.size(); i++)
		zones[i].resize(rows);

	for (int i = 0; i < cols; i++)
	{
		const int start_width = floor((double)i * zone_width);
		const int stop_width = floor(double(i + 1) * zone_width);
		for (int j = 0; j < rows; j++)
		{
			const int start_height = floor((double)j * zone_height);
			const int stop_height = floor(double(j + 1) * zone_height);
			const frequencies freq = GetGScFrequency(start_width, start_height, stop_width, stop_height);
			const int area = (stop_width - start_width) * (stop_height - start_height);

			zones[i][j] = std::unique_ptr<zonestat>(new zonestat(Point(i, j), freq));

			//Shannon entropy with 32 different values:
			constexpr int group_size = 8;
			int contract_freq[256 / group_size] = { 0 };
			double prob[256 / group_size] = { 0 };

			int most_populated = 0;
			for (int k = 0; k < 256; k++)
			{
				contract_freq[k / group_size] += freq[k];
				if (contract_freq[k / group_size] >= contract_freq[most_populated])
					most_populated = k / group_size;
			}
			double entropy = 0;
			for (int k = 0; k < (256 / group_size); k++)
			{
				if (k == most_populated) continue;
				prob[k] = (double)contract_freq[k] / area;
				if (prob[k] > 0)
					entropy += prob[k] * log2(prob[k]);
			}

			zones[i][j]->entropy = entropy * -1;
			if (zones[i][j]->entropy > zonestat::max_entropy)
				zonestat::max_entropy = zones[i][j]->entropy;
		}
	}

	/* Finding reference zone and mapping entropy values */
	reference_zone = zones[0][0].get();
	std::map<ubyte, unsigned short> entropy_freq;

	for (int i = 0; i < cols; i++)
	{
		for (int j = 0; j < rows; j++)
		{
			if (zones[i][j]->peak > reference_zone->peak)
				reference_zone = zones[i][j].get();

			/* Mapping to 0.0 - 1.0 */
			zones[i][j]->relative_entropy = zones[i][j]->entropy / zonestat::max_entropy;
			entropy_freq[zones[i][j]->relative_entropy]++;
			if (entropy_freq[zones[i][j]->relative_entropy] > zonestat::most_common_relative_entropy)
				zonestat::most_common_relative_entropy = zones[i][j]->relative_entropy;

			/* Mapping to 0 - 255 */
			ubyte entropy_shade = floor(zones[i][j]->relative_entropy * 255);
			zones[i][j]->entropy_byte = entropy_shade;
			entropy_freq[entropy_shade]++;
			if (entropy_freq[entropy_shade] > zonestat::most_common_entropy_byte)
				zonestat::most_common_entropy_byte = entropy_shade;
		}
	}

	return amount;
}

void Document::ScaleLighting(bool scale_rgb)
{
	/* Scaling each pixel: */
	double** pixel_scale_factor = new double* [width];
	for (size_t i = 0; i < width; i++)
	{
		pixel_scale_factor[i] = new double[height];
		for (size_t j = 0; j < height; j++)
		{
			pixel_scale_factor[i][j] = 1.0;
		}
	}

	auto triangles = Triangulation(ZonesToTriangulate());

	for (int i = 0; i < triangles.size(); i++)
	{
		//Triangle:
		Triangle triangle = triangles[i];
		Plane* plane = &triangle.plane;
		Point* A = &triangle.A;
		Point* B = &triangle.B;
		Point* C = &triangle.C;

		//Bounding rectangle of the triangle:
		Point TL = Point(*A);
		//x:
		if (B->x < TL.x) TL.x = B->x;
		if (C->x < TL.x) TL.x = C->x;
		//y:
		if (B->y < TL.y) TL.y = B->y;
		if (C->y < TL.y) TL.y = C->y;
		Point BR = Point(*A);
		//x:
		if (B->x > BR.x) BR.x = B->x;
		if (C->x > BR.x) BR.x = C->x;
		//y:
		if (B->y > BR.y) BR.y = B->y;
		if (C->y > BR.y) BR.y = C->y;

		for (size_t x = TL.x; x <= BR.x; x++)
		{
			for (size_t y = TL.y; y <= BR.y; y++)
			{
				double* coords = CartesianToBarycentric(x, y, *A, *B, *C);

				if ((0 <= coords[0] && coords[0] <= 1) &&
					(0 <= coords[1] && coords[1] <= 1) &&
					(0 <= coords[2] && coords[2] <= 1))	//The coordinates are inside of the triangle
				{
					/* Planar scaling: */
					pixel_scale_factor[x][y] = plane->Interpolate(x, y);
				}
				delete[] coords;
			}
		}
	}

	for (size_t i = 0; i < width; i++)
	{
		for (size_t j = 0; j < height; j++)
		{
			if (pixel_scale_factor[i][j] <= 1.0)
			{
				int left_ind = i;
				while ((left_ind >= 0) && (pixel_scale_factor[left_ind][j] <= 1.0))
					left_ind--;
				int right_ind = i;
				while ((right_ind < width) && (pixel_scale_factor[right_ind][j] <= 1.0))
					right_ind++;

				int up_ind = j;
				while ((up_ind >= 0) && (pixel_scale_factor[i][up_ind] <= 1.0))
					up_ind--;
				int down_ind = j;
				while ((down_ind < height) && (pixel_scale_factor[i][down_ind] <= 1.0))
					down_ind++;

				bool x_wise = (right_ind - left_ind) <= (down_ind - up_ind);
				if (left_ind < 0 || right_ind >= width) x_wise = false;
				else if (up_ind < 0 || down_ind >= height) x_wise = true;

				if (x_wise)
				{
					for (size_t k = left_ind + 1; k < right_ind; k++)
					{
						double tmp = (k - left_ind) * pixel_scale_factor[right_ind][j] + (right_ind - k) * pixel_scale_factor[left_ind][j];
						tmp /= (right_ind - left_ind);
						pixel_scale_factor[k][j] = tmp;
					}
				}
				else
				{
					for (size_t k = up_ind + 1; k < down_ind; k++)
					{
						double tmp = (k - up_ind) * pixel_scale_factor[i][down_ind] + (down_ind - k) * pixel_scale_factor[i][up_ind];
						tmp /= (down_ind - up_ind);
						pixel_scale_factor[i][k] = tmp;
					}
				}
			}

			int new_val = gscpixels[i][j] * pixel_scale_factor[i][j];
			if (new_val > GScPixel::WHITE)
				new_val = GScPixel::WHITE;
			if (new_val < GScPixel::BLACK)
				//new_val = GScPixel::BLACK;
				new_val = gscpixels[i][j];
			gscpixels[i][j] = new_val;

			if (scale_rgb)
			{
				int r = rgbpixels[i][j].red * pixel_scale_factor[i][j];
				if (r > 0xFF) r = 0xFF;
				int g = rgbpixels[i][j].green * pixel_scale_factor[i][j];
				if (g > 0xFF) g = 0xFF;
				int b = rgbpixels[i][j].blue * pixel_scale_factor[i][j];
				if (b > 0xFF) b = 0xFF;

				rgbpixels[i][j].red = r;
				rgbpixels[i][j].green = g;
				rgbpixels[i][j].blue = b;
			}
		}
	}

	for (size_t i = 0; i < width; i++)
		delete[] pixel_scale_factor[i];
	delete[] pixel_scale_factor;
}

std::vector<zonestat*> Document::ZonesToTriangulate()
{
	const zonestat global_stat = zonestat(GetGScFrequency());

	std::vector<zonestat*> triangulation_zones;
	//Adding every zone that has the most common entropy value and is not too dark to triangulation_zones
	for (int i = 0; i < zonestat::zonestat::cols; i++)
	{
		const int start_width = floor((double)i * zonestat::width);
		const int stop_width = floor(double(i + 1) * zonestat::width);

		for (int j = 0; j < zonestat::zonestat::rows; j++)
		{
			if ((zones[i][j]->entropy_byte <= zonestat::most_common_entropy_byte) &&
				(zones[i][j]->whitest >= 0.9 * (global_stat.peak - (GScPixel::WHITE - global_stat.peak)))) //If the whitest pixel in the zone is too dark then it's brightness would scale down, not up.
			{
				const int start_height = floor((double)j * zonestat::height);
				const int stop_height = floor(double(j + 1) * zonestat::height);

				const double mult = (double)reference_zone->peak / zones[i][j]->peak;

				zones[i][j]->scale_factor = mult;
				triangulation_zones.push_back(zones[i][j].get());
			}
		}
	}

	//Finding the coordinates for the outer most points in triangulation_zones
	int upper = zonestat::rows - 1, lower = 0, leftmost = zonestat::cols - 1, rightmost = 0;
	size_t closest_TL = 0, closest_BL = 0, closest_TR = 0, closest_BR = 0;
	for (size_t i = 0; i < triangulation_zones.size(); i++)
	{
		if (triangulation_zones[i]->point->x <= leftmost)
		{
			leftmost = triangulation_zones[i]->point->x;
			if (Point::DistanceSq(*triangulation_zones[i]->point, *zones[0][0]->point) <
				Point::DistanceSq(*triangulation_zones[i]->point, *triangulation_zones[closest_TL]->point))
				closest_TL = i;
			if (Point::DistanceSq(*triangulation_zones[i]->point, *zones[0][zonestat::rows - 1]->point) <
				Point::DistanceSq(*triangulation_zones[i]->point, *triangulation_zones[closest_BL]->point))
				closest_BL = i;
		}
		if (triangulation_zones[i]->point->x >= rightmost)
		{
			rightmost = triangulation_zones[i]->point->x;
			if (Point::DistanceSq(*triangulation_zones[i]->point, *zones[zonestat::cols - 1][0]->point) <
				Point::DistanceSq(*triangulation_zones[i]->point, *triangulation_zones[closest_TR]->point))
				closest_TR = i;
			if (Point::DistanceSq(*triangulation_zones[i]->point, *zones[zonestat::cols - 1][zonestat::rows - 1]->point) <
				Point::DistanceSq(*triangulation_zones[i]->point, *triangulation_zones[closest_BR]->point))
				closest_BR = i;
		}
		if (triangulation_zones[i]->point->y <= upper)
		{
			upper = triangulation_zones[i]->point->y;
			if (Point::DistanceSq(*triangulation_zones[i]->point, *zones[0][0]->point) <
				Point::DistanceSq(*triangulation_zones[i]->point, *triangulation_zones[closest_TL]->point))
				closest_TL = i;
			if (Point::DistanceSq(*triangulation_zones[i]->point, *zones[zonestat::cols - 1][0]->point) <
				Point::DistanceSq(*triangulation_zones[i]->point, *triangulation_zones[closest_TR]->point))
				closest_TR = i;
		}
		if (triangulation_zones[i]->point->y >= lower)
		{
			lower = triangulation_zones[i]->point->y;
			if (Point::DistanceSq(*triangulation_zones[i]->point, *zones[0][zonestat::rows - 1]->point) <
				Point::DistanceSq(*triangulation_zones[i]->point, *triangulation_zones[closest_BL]->point))
				closest_BL = i;
			if (Point::DistanceSq(*triangulation_zones[i]->point, *zones[zonestat::cols - 1][zonestat::rows - 1]->point) <
				Point::DistanceSq(*triangulation_zones[i]->point, *triangulation_zones[closest_BR]->point))
				closest_BR = i;
		}
	}
	//Finding a convex hull-like line that contains the outer most points of the zones (but instead of being minimal, it's maximal)
	std::vector<int> hull_points_ind;
	Point TL = *triangulation_zones[0]->point, BL = *triangulation_zones[0]->point,
		BR = *triangulation_zones[0]->point, TR = *triangulation_zones[0]->point;
	for (size_t i = 0; i < triangulation_zones.size(); i++)
	{
		Point curr = *triangulation_zones[i]->point;

		if (curr.x == leftmost)
		{
			if (Point::DistanceSq(curr, Point(0, 0)) < Point::DistanceSq(Point(0, 0), TL))
				TL = curr;
			if (Point::DistanceSq(curr, Point(0, zonestat::rows - 1)) < Point::DistanceSq(Point(0, zonestat::rows - 1), BL))
				BL = curr;
			hull_points_ind.push_back(i);
		}
		else if (curr.y == lower)
		{
			if (Point::DistanceSq(curr, Point(0, zonestat::rows - 1)) < Point::DistanceSq(Point(0, zonestat::rows - 1), BL))
				BL = curr;
			if (Point::DistanceSq(curr, Point(zonestat::cols - 1, zonestat::rows - 1)) < Point::DistanceSq(Point(zonestat::cols - 1, zonestat::rows - 1), BR))
				BR = curr;
			hull_points_ind.push_back(i);
		}
		else if (curr.x == rightmost)
		{
			if (Point::DistanceSq(curr, Point(zonestat::cols - 1, 0)) < Point::DistanceSq(Point(zonestat::cols - 1, 0), TR))
				TR = curr;
			if (Point::DistanceSq(curr, Point(zonestat::cols - 1, zonestat::rows - 1)) < Point::DistanceSq(Point(zonestat::cols - 1, zonestat::rows - 1), BR))
				BR = curr;
			hull_points_ind.push_back(i);
		}
		else if (curr.y == upper)
		{
			if (Point::DistanceSq(curr, Point(0, 0)) < Point::DistanceSq(Point(0, 0), TL))
				TL = curr;
			if (Point::DistanceSq(curr, Point(zonestat::cols - 1, 0)) < Point::DistanceSq(Point(zonestat::cols - 1, 0), TR))
				TR = curr;
			hull_points_ind.push_back(i);
		}
	}
	//Finding the four outer corners of the found zones: (based on distance of the overall corners)
	for (size_t i = 0; i < hull_points_ind.size(); i++)
	{
		Point* curr = triangulation_zones[hull_points_ind[i]]->point.get();

		if (curr->equals(TL))
		{
			if (Point::DistanceSq(Point(0, 0), TL) > 0)
			{
				zones[0][0]->scale_factor = curr->zone->scale_factor;
				triangulation_zones.push_back(zones[0][0].get());
			}
			continue;
		}
		else if (curr->equals(TR))
		{
			if (Point::DistanceSq(Point(zonestat::cols - 1, 0), TR) > 0)
			{
				zones[zonestat::cols - 1][0]->scale_factor = curr->zone->scale_factor;
				triangulation_zones.push_back(zones[zonestat::cols - 1][0].get());
			}
			continue;
		}
		else if (curr->equals(BL))
		{
			if (Point::DistanceSq(Point(0, zonestat::rows - 1), BL) > 0)
			{
				zones[0][zonestat::rows - 1]->scale_factor = curr->zone->scale_factor;
				triangulation_zones.push_back(zones[0][zonestat::rows - 1].get());
			}
			continue;
		}
		else if (curr->equals(BR))
		{
			if (Point::DistanceSq(Point(zonestat::cols - 1, zonestat::rows - 1), BR) > 0)
			{
				zones[zonestat::cols - 1][zonestat::rows - 1]->scale_factor = curr->zone->scale_factor;
				triangulation_zones.push_back(zones[zonestat::cols - 1][zonestat::rows - 1].get());
			}
			continue;
		}

		if (curr->x == leftmost)
		{
			zones[0][curr->y]->scale_factor = curr->zone->scale_factor;
			triangulation_zones.push_back(zones[0][curr->y].get());
		}
		else if (curr->y == lower)
		{
			zones[curr->x][zonestat::rows - 1]->scale_factor = curr->zone->scale_factor;
			triangulation_zones.push_back(zones[curr->x][zonestat::rows - 1].get());
		}
		else if (curr->x == rightmost)
		{
			zones[zonestat::cols - 1][curr->y]->scale_factor = curr->zone->scale_factor;
			triangulation_zones.push_back(zones[zonestat::cols - 1][curr->y].get());
		}
		else if (curr->y == upper)
		{
			zones[curr->x][0]->scale_factor = curr->zone->scale_factor;
			triangulation_zones.push_back(zones[curr->x][0].get());
		}
	}

	return triangulation_zones;
}

std::vector<Triangle> Document::Triangulation(std::vector<zonestat*> zones)
{
	/* Meshing: */

	std::vector<CDT::V2d<double>> triangulation_points;
	//Converting custom Point struct to the one used by the library:
	for (size_t i = 0; i < zones.size(); i++)
	{
		Point* P = zones[i]->point.get();
		triangulation_points.push_back(CDT::V2d<double>::make(P->x, P->y));
	}

	CDT::RemoveDuplicates(triangulation_points);
	CDT::Triangulation<double> triangulation;
	triangulation.insertVertices(triangulation_points);
	triangulation.eraseSuperTriangle();

	//Converting zone coordinates to pixel coordinates:
	std::vector<Point> vertex_pixels;
	for (auto point = triangulation_points.begin(); point < triangulation_points.end(); point++)
	{
		int x = (int)(*point).x;
		int y = (int)(*point).y;
		//Calculating central pixel of the zone:
		Point tmp = Point(x * zonestat::width + zonestat::width / 2, y * zonestat::height + zonestat::height / 2, this->zones[x][y].get());
		//If the zone is on a side the corresponding coordinate must be too.
		//(otherwise a border will be left on the image, if corener and edge zones are also converted to central pixels)
		if ((*point).x == 0)
			tmp.x = 0;
		else if ((*point).x == zonestat::cols - 1)
			tmp.x = width - 1;
		if ((*point).y == 0)
			tmp.y = 0;
		else if ((*point).y == zonestat::rows - 1)
			tmp.y = height - 1;
		vertex_pixels.push_back(tmp);
	}
	//Creating the triangles:
	std::vector<Triangle> triangles;
	for (size_t i = 0; i < triangulation.triangles.size(); i++)
	{
		auto current = triangulation.triangles[i].vertices;
		triangles.push_back(Triangle(vertex_pixels[current[0]], vertex_pixels[current[1]], vertex_pixels[current[2]]));
	}

	return triangles;
}

void Document::AdjustContrast(const char adjustment)
{
	/* Adjusting contrast: */
	double factor = (259 * double(adjustment + 255)) / (255 * (259 - adjustment));
	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			int newval = factor * (gscpixels[i][j] - adjustment) + adjustment;
			if (newval < 0) newval = 0;
			if (newval > 255) newval = 255;
			gscpixels[i][j] = newval;
		}
	}
}

void Document::DeleteBackground()
{
	if (gscpixels.empty()) RGBtoGSc();
	if (zones.empty()) InitZones();
	
	/* Background removal: */
	constexpr double BG = 0.6;
	constexpr double TEXT = 0.7;

	unsigned int BGfreq[256]{ 0 };
	unsigned long int BGarea = 0;
	unsigned int TEXTfreq[256]{ 0 };
	unsigned long int TEXTarea = 0;

	for (int i = 0; i < zonestat::cols; i++)
	{
		const int start_width = floor((double)i * zonestat::width);
		const int stop_width = floor(double(i + 1) * zonestat::width);
		for (int j = 0; j < zonestat::rows; j++)
		{
			const int start_height = floor((double)j * zonestat::height);
			const int stop_height = floor(double(j + 1) * zonestat::height);

			const frequencies freq = GetGScFrequency(start_width, start_height, stop_width, stop_height);
			const zonestat zone = zonestat(freq);
			const unsigned int area = (stop_width - start_width) * (stop_height - start_height);

			if ((zones[i][j]->relative_entropy < BG) && (zone.peak >= 30))
			{
				for (int k = 0; k < 256; k++)
				{
					BGfreq[k] += freq[k];
					BGarea += area;
				}
			}
			else if (zones[i][j]->relative_entropy > TEXT)
			{
				for (int k = 0; k < 256; k++)
				{
					TEXTfreq[k] += freq[k];
					TEXTarea += area;
				}
			}
		}
	}

	const double area_norming = (double)TEXTarea / BGarea;
	double prob[256]{ 0.0 };

	for (size_t i = 0; i < 20; i++)
		prob[i] = -1;
	for (size_t i = 20; i < 230; i++)
	{
		prob[i] = (((double)BGfreq[i] * area_norming) - (double)TEXTfreq[i]);
	}
	for (size_t i = 20; i < 230 - 1; i++)
	{
		if (prob[i - 1] < 0 && prob[i] > 0 && prob[i + 1] < 0)
			prob[i] = -1;
	}
	for (size_t i = 230; i <= GScPixel::WHITE; i++)
		prob[i] = 1;

	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			ubyte& shade = gscpixels[i][j];

			if (prob[shade] > 0)
				shade = GScPixel::WHITE;
		}
	}

	is_background_removed = true;
}

bool Document::Read()
{
	ilGenImages(1, &image_name);
	ilBindImage(image_name);
	bool success = ilLoadImage(image_path.c_str());
	ILenum Error = ilGetError();
	if (success)
	{
		width = ilGetInteger(IL_IMAGE_WIDTH);
		height = ilGetInteger(IL_IMAGE_HEIGHT);

		ILDataToRGB();
	}
	return success;
}

bool Document::WriteRGB(std::string file_to_create, int file_type) const
{
	ILuint new_name;
	ilGenImages(1, &new_name);
	ilBindImage(new_name);

	ILubyte* new_data = RGBToILData(true);
	bool success = false;
	success = ilTexImage(width, height, 0, 4, IL_RGBA, IL_UNSIGNED_BYTE, new_data);
	ilRegisterOrigin(IL_ORIGIN_UPPER_LEFT);
	//ILenum Error = ilGetError();
	delete[] new_data;

	if (success)
	{
		if (file_type == -1)
			success = ilSaveImage(file_to_create.c_str());
		else
			success = ilSave(file_type, file_to_create.c_str());
	}

	//ILenum Error = ilGetError();
	ilDeleteImage(new_name);

	return success;
}

bool Document::WriteGSc(std::string file_to_create, int file_type)
{
	if (gscpixels.empty()) RGBtoGSc();

	ILuint new_name;
	ilGenImages(1, &new_name);
	ilBindImage(new_name);

	ILubyte* new_data = GScToILData();
	bool success = false;
	success = ilTexImage(width, height, 0, 1, IL_LUMINANCE, IL_UNSIGNED_BYTE, new_data);
	ilRegisterOrigin(IL_ORIGIN_UPPER_LEFT);
	//ILenum Error = ilGetError();
	delete[] new_data;

	if (success)
	{
		if (file_type == -1)
			success = ilSaveImage(file_to_create.c_str());
		else
			success = ilSave(file_type, file_to_create.c_str());
	}

	//ILenum Error = ilGetError();
	ilDeleteImage(new_name);

	return success;
}