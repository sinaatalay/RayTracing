#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

//==========CLASSES=====================================
class vec {  //Class for vectors and RGB colors
public:
	//Attributes
	double x = 0, y = 0, z = 0;

	//Constructors
	vec() {}
	vec(double a, double b, double c) : x(a), y(b), z(c) {}

	// Basic operations
	vec operator-(const vec& a) const {  //Vector difference
		vec result;
		result.x = x - a.x;
		result.y = y - a.y;
		result.z = z - a.z;
		return result;
	}
	vec operator+(const vec& a) const {  //Vector addition
		vec result;
		result.x = x + a.x;
		result.y = y + a.y;
		result.z = z + a.z;
		return result;
	}
	vec operator*(double a) const {  //Scalar multiplication
		vec result;
		result.x = x * a;
		result.y = y * a;
		result.z = z * a;
		return result;
	}
	double operator*(const vec& a) const {  //Dot product
		double result;
		result = x * a.x + y * a.y + z * a.z;
		return result;
	}
	vec operator&(const vec& b) const {  //Cross product
		vec result;
		result.x = y * b.z - z * b.y;
		result.y = z * b.x - x * b.z;
		result.z = x * b.y - y * b.x;
		return result;
	}
	vec operator/(double a) const {  //Scalar division
		vec result;
		result.x = x / a;
		result.y = y / a;
		result.z = z / a;
		return result;
	}
	double length() const {  //Length of the vector
		return std::sqrt(x * x + y * y + z * z);
	}
	vec unitVector() const {  //Unit direction of the vector
		vec result;
		double length = (*this).length();
		result.x = x / length;
		result.y = y / length;
		result.z = z / length;
		return result;
	}
};
using color = vec;                       //To not confuse vectors with colors
vec operator*(double a, const vec& b) {  //To generalize scalar multiplication
	vec result;
	result = b * a;
	return result;
}
vec operator/(double a, const vec& b) {  //To generalize scalar division
	vec result;
	result = b / a;
	return result;
}

class sphere {  //Class for spheres
public:
	//Attributes
	double radius;
	vec center, color, coeffs;

	//Constructor
	sphere(const vec& a, double b, const vec& c, double kA, double kD, double kS) {
		center = a;
		radius = b;
		color = c;
		coeffs.x = kA;  //Ambient term
		coeffs.y = kD;  //Diffuse term
		coeffs.z = kS;  //Specular term
	}
};

class plane {  //Class for planes
public:
	//Attributes
	vec normal, point, color, coeffs;

	//Constructor
	plane(const vec& a, const vec& b, const vec& c, double kA, double kD, double kS) {
		normal = a.unitVector();
		point = b;
		color = c;
		coeffs.x = kA;  //Ambient term
		coeffs.y = kD;  //Diffuse term
		coeffs.z = kS;  //Specular term
	}
};

// Triangles are not used for this scene:
class triangle {
public:
	//Attributes
	vec normal, cm, ver1, ver2, ver3, color, coeffs;

	//Constructor
	triangle(const vec& b1, const vec& b2, const vec& b3, const vec& c, double kA, double kD, double kS) {
		vec first = b2 - b1;
		vec second = b3 - b1;
		cm = (b1 + b2 + b3) / 3;
		normal = (first & second).unitVector();
		ver1 = b1;
		ver2 = b2;
		ver3 = b3;
		color = c;
		coeffs.x = kA;  //Ambient term
		coeffs.y = kD;  //Diffuse term
		coeffs.z = kS;  //Specular term
	}

	bool isInside(const vec& point) const {
		if (((ver2 - ver1) & (point - ver1)) * ((ver2 - ver1) & (cm - ver1)) < 0.0000001) {
			return 0;
		}
		else if (((ver3 - ver1) & (point - ver1)) * ((ver3 - ver1) & (cm - ver1)) < 0.0000001) {
			return 0;
		}
		else if (((ver2 - ver3) & (point - ver3)) * ((ver2 - ver3) & (cm - ver3)) < 0.0000001) {
			return 0;
		}
		else {
			return 1;
		}
	}
};

class ray {  //Class for rays
public:
	//Attributes
	vec start, direction;

	//Constructors
	ray() {}
	ray(const vec& a, const vec& b) {
		start = a;
		direction = b;
	}

	//Basic operations
	vec where(double t) const {  //Calculate where = start + direction*t
		vec result;
		result = start + t * direction;
		return result;
	}
};

class lightSource {
public:
	//Attributes
	vec position;
	double I_a, I_source;

	//Constructor
	lightSource(const vec& a, double I_ambient, double I_lightSource) {
		position = a;
		I_a = I_ambient;
		I_source = I_lightSource;
	}
};
//======================================================

//==========FUNCTIONS===================================
bool isShadow(const vec& intersectionPoint, const lightSource& light, const std::vector<sphere>& spheres, const std::vector<plane>& planes) {
	//Create a ray from intersection point to light source:
	ray lightRay(intersectionPoint, light.position - intersectionPoint);

	//First check intersections with spheres
	std::size_t count = spheres.size();
	for (std::size_t i = 0; i < count; i++) {  //Try for all spheres
		double a = lightRay.direction * lightRay.direction;
		double b = 2 * (lightRay.direction * (lightRay.start - spheres[i].center));
		double c = (lightRay.start - spheres[i].center) * (lightRay.start - spheres[i].center) - spheres[i].radius * spheres[i].radius;
		double discriminant = (b * b - 4 * a * c);
		if (discriminant >= 0) {
			double sqrtDiscriminant = std::sqrt(discriminant);
			double firstSolution = (-b - sqrtDiscriminant) / (2 * a);
			if (firstSolution > 0.001) {  //Solutions with negative t doesn't help us
				return true;
			}
			else {
				double secondSolution = (-b + sqrtDiscriminant) / (2 * a);
				if (secondSolution > 0.001) {  //Solutions with negative t doesn't help us
					return true;
				}
			}
		}
	}

	//NOTE: For this scene, planes have no effect on shadows.
	//Then check intersections with planes
	// count = planes.size();
	// for (std::size_t i = 0; i < count; i++) {  //Try for all planes
	//     double a = planes[i].normal * (planes[i].point - lightRay.start);
	//     double b = lightRay.direction * planes[i].normal;
	//     double c = a / b;
	//     if (c > 0.001) {
	//         return true;
	//     }
	// }

	//NOTE: For this scene, triangles have no effect on shadows.
	// count = triangles.size();
	// for (std::size_t i = 0; i < count; i++) {  //Try for all triangles
	//	   double a = triangles[i].normal * (triangles[i].ver1 - lightRay.start);
	//	   double b = lightRay.direction * triangles[i].normal;
	//	   double c = a / b;
	//	   vec Tpoint = lightRay.where(c);
	//	   if (c > 0.001 && triangles[i].isInside(Tpoint)) {
	//	   	return true;
	//	   }
	// }

	//If there are no intersection
	return false;
}
color calculateColor(const std::vector<vec>& intersection, const lightSource& light, const std::vector<sphere>& spheres, const std::vector<plane>& planes, const vec& currentDir) {
	color ambient, diffusive, specular, defaultColor = intersection[1], white(255, 255, 255), result;
	double kA = intersection[3].x, kD = intersection[3].y, kS = intersection[3].z;
	vec relativePos = light.position - intersection[0];
	vec l = relativePos.unitVector();
	ambient = kA * light.I_a * defaultColor;
	bool shadow = isShadow(intersection[0], light, spheres, planes);
	if (!shadow) {
		double distance = relativePos.length();
		diffusive = kD * (light.I_source / (distance * distance)) * std::max(0.0, l * intersection[2]) * defaultColor;
		vec reflection = (-1 * l) - 2 * (-1 * l * intersection[2]) * intersection[2];
		specular = kS * (light.I_source / (distance * distance)) * pow(std::max(0.0, reflection * (-1 * currentDir).unitVector()), 20) * white;
	}
	result = ambient + diffusive + specular;
	result.x = std::min(255.0, result.x);
	result.y = std::min(255.0, result.y);
	result.z = std::min(255.0, result.z);
	return result;
}
std::vector<vec> findClosestIntersection(const ray& ray, const std::vector<sphere>& spheres, const std::vector<plane>& planes) {
	std::vector<vec> result;
	int closestSphere = 0, closestPlane = 0, closestTriangle = 0;
	double smallestSolutionSphere = 50000, smallestSolutionPlane = 50000, smallestSolutionTriangle = 50000;

	//First check intersections with spheres
	//Solution of at^2+bt+c=0 type of equation
	std::size_t count = spheres.size();
	std::vector<double> solutions, index;
	for (std::size_t i = 0; i < count; i++) {  //Try for all spheres
		double a = ray.direction * ray.direction;
		double b = 2 * (ray.direction * (ray.start - spheres[i].center));
		double c = (ray.start - spheres[i].center) * (ray.start - spheres[i].center) - spheres[i].radius * spheres[i].radius;
		double discriminant = (b * b - 4 * a * c);
		if (discriminant >= 0) {
			double sqrtDiscriminant = std::sqrt(discriminant);
			double firstSolution = (-b - sqrtDiscriminant) / (2 * a);
			if (firstSolution > 0.001 && firstSolution < 50000) {  //Solutions with negative t doesn't help us
				index.push_back(i);
				solutions.push_back(firstSolution);
			}
			else {
				double secondSolution = (-b + sqrtDiscriminant) / (2 * a);
				if (secondSolution > 0.001 && secondSolution < 50000) {  //Solutions with negative t doesn't help us
					index.push_back(i);
					solutions.push_back(secondSolution);
				}
			}
		}
	}
	//Find the smallest 't' to find closest sphere
	if (solutions.size() > 0) {
		double temp = solutions[0];
		int j = 0;
		for (std::size_t i = 1; i < solutions.size(); i++) {
			if (solutions[i] < temp) {
				temp = solutions[i];
				j = i;
			}
		}
		//Store closest sphere and solution t
		closestSphere = index[j];
		smallestSolutionSphere = solutions[j];
	}

	//======================================================

	//Then check intersections with planes
	count = planes.size();
	std::vector<double> solutionsPlanes, indexPlanes;
	for (std::size_t i = 0; i < count; i++) {  //Try for all planes
		double a = planes[i].normal * (planes[i].point - ray.start);
		double b = ray.direction * planes[i].normal;
		double c = a / b;
		if (c > 0.001 && c < 50000) {
			indexPlanes.push_back(i);
			solutionsPlanes.push_back(c);
		}
	}
	//Find the smallest 't' to find closest plane
	if (solutionsPlanes.size() > 0) {
		double temp = solutionsPlanes[0];
		int j = 0;
		for (std::size_t i = 1; i < solutionsPlanes.size(); i++) {
			if (solutionsPlanes[i] < temp) {
				temp = solutionsPlanes[i];
				j = i;
			}
		}
		//Store closest plane and solution t
		closestPlane = indexPlanes[j];
		smallestSolutionPlane = solutionsPlanes[j];
	}

	//======================================================

	//Compare the closest plane and sphere
	if (smallestSolutionPlane < smallestSolutionSphere) {
		result.push_back(ray.where(smallestSolutionPlane));  //Surface intersection point
		result.push_back(planes[closestPlane].color);        //Surface default color
		result.push_back(planes[closestPlane].normal);       //Surface normal unit vector.
		result.push_back(planes[closestPlane].coeffs);       //Ambient and diffusive coefficients.
		result.push_back(planes[closestPlane].coeffs);       //Dummy to see if it's plane or not.
	}
	else if (smallestSolutionPlane > smallestSolutionSphere) {
		vec point = ray.where(smallestSolutionSphere);
		result.push_back(point);                                                 //Surface intersection point
		result.push_back(spheres[closestSphere].color);                          //Surface default color
		result.push_back((point - spheres[closestSphere].center).unitVector());  //Surface normal unit vector.
		result.push_back(spheres[closestSphere].coeffs);                         //Ambient and diffusive coefficients.
	}
	return result;
}
color findColor(const ray& firstRay, const std::vector<sphere>& spheres, const std::vector<plane>& planes, int maxDepth, const lightSource& light) {
	color colorResult, noIntersection(255, 255, 255);
	int depth = 0;
	ray currentRay = firstRay;
	while (depth < maxDepth) {
		std::vector<vec> intersection = findClosestIntersection(currentRay, spheres, planes);
		if (intersection.size() == 5) { //If plane, no reflection
			maxDepth = 1;
		}
		if (intersection.size() > 0) {
			depth = depth + 1;
			colorResult = colorResult + calculateColor(intersection, light, spheres, planes, currentRay.direction);
			currentRay.direction = currentRay.direction - 2 * (currentRay.direction * intersection[2]) * intersection[2];
			currentRay.start = intersection[0];
		}
		else {
			break;
		}
	}
	if (depth != 0) {
		colorResult = colorResult / depth;
	}
	else {
		return noIntersection;
	}
	return colorResult;
}
//======================================================

//==========BMP FILE====================================
std::ofstream createBmpFile(int H, int W) {
	//To create BMP file
	//BMP code from https://youtu.be/vqT5j38bWGg
	std::ofstream image;
	image.open("output.bmp", std::ios::binary | std::ios::out);
	const int paddingAmount = ((4 - (W * 3) % 4) % 4);
	const int fileHeaderSize = 14;
	const int informationHeaderSize = 40;
	const int fileSize = fileHeaderSize + informationHeaderSize + W * H * 3 + paddingAmount * H;

	unsigned char fileHeader[fileHeaderSize];
	for (int i = 0; i <= 13; i++) {
		fileHeader[i] = 0;
	}
	fileHeader[0] = 'B';
	fileHeader[1] = 'M';
	fileHeader[2] = fileSize;
	fileHeader[3] = fileSize >> 8;
	fileHeader[4] = fileSize >> 16;
	fileHeader[5] = fileSize >> 24;
	fileHeader[10] = fileHeaderSize + informationHeaderSize;

	unsigned char informationHeader[informationHeaderSize];
	for (int i = 0; i <= 39; i++) {
		informationHeader[i] = 0;
	}
	informationHeader[0] = informationHeaderSize;
	informationHeader[4] = W;
	informationHeader[5] = W >> 8;
	informationHeader[6] = W >> 16;
	informationHeader[7] = W >> 24;
	informationHeader[8] = H;
	informationHeader[9] = H >> 8;
	informationHeader[10] = H >> 16;
	informationHeader[11] = H >> 24;
	informationHeader[12] = 1;
	informationHeader[14] = 24;

	image.write(reinterpret_cast<char*>(fileHeader), fileHeaderSize);
	image.write(reinterpret_cast<char*>(informationHeader), informationHeaderSize);

	return image;
}
void writeColor(color pixelColor, std::ofstream& image) {  //Write RGB color of a pixel to the BMP file
	unsigned char r = static_cast<unsigned char>(pixelColor.x);
	unsigned char g = static_cast<unsigned char>(pixelColor.y);
	unsigned char b = static_cast<unsigned char>(pixelColor.z);
	unsigned char color[] = { b, g, r };
	image.write(reinterpret_cast<char*>(color), 3);
}
//======================================================

int main() {
	//In this code, ambient, diffuse, and specular coefficients are taken 0.3, 0.4, and 0.3 respectively; for all objects.
	double kA = 0.3, kD = 0.4, kS = 0.3;
	int maxDepth = 3;

	std::vector<sphere> spheres;
	vec sphereV1(-36, -28, 220);
	vec sphereV2(-55, 23, 230);
	vec sphereV3(-107, -39, 157);
	vec sphereV4(13, -11, 235);
	vec sphereV5(-151, -27, 220);
	vec sphereV6(-81, -32, 160);
	vec sphereV7(7, -36, 164);
	vec sphereV8(-46, -37, 179);
	vec sphereV9(-5, -11, 278);
	vec sphereV10(6, -32, 199);
	vec sphereV11(-24, -28, 197);
	vec sphereV12(34, -20, 232);
	vec sphereV13(-120, -2, 333);
	vec sphereV14(-126, -27, 199);
	vec sphereV15(-99, -11, 291);
	vec sphereV16(53, -5, 288);
	vec sphereV17(-120, -20, 260);
	vec sphereV18(-145, -15, 263);
	vec sphereV19(-153, 96, 426);
	vec sphereV20(5, 19, 330);
	color sphereC1(139, 0, 0);
	color sphereC2(255, 0, 0);
	color sphereC3(255, 69, 0);
	color sphereC4(255, 215, 0);
	color sphereC5(128, 128, 0);
	color sphereC6(255, 255, 0);
	color sphereC7(124, 252, 0);
	color sphereC8(0, 100, 0);
	color sphereC9(0, 128, 0);
	color sphereC10(0, 255, 0);
	color sphereC11(0, 255, 255);
	color sphereC12(0, 191, 255);
	color sphereC13(0, 0, 128);
	color sphereC14(65, 105, 225);
	color sphereC15(139, 9, 139);
	color sphereC16(128, 0, 255);
	color sphereC17(230, 230, 250);
	color sphereC18(112, 128, 144);
	color sphereC19(176, 196, 222);
	color sphereC20(192, 192, 192);
	double radius1 = 10;
	double radius2 = 35;
	double radius3 = 10;
	double radius4 = 10;
	double radius5 = 10;
	double radius6 = 14;
	double radius7 = 10;
	double radius8 = 10;
	double radius9 = 10;
	double radius10 = 15;
	double radius11 = 13;
	double radius12 = 11;
	double radius13 = 18;
	double radius14 = 12;
	double radius15 = 8;
	double radius16 = 10;
	double radius17 = 10;
	double radius18 = 10;
	double radius19 = 81;
	double radius20 = 35;
	sphere sphere1(sphereV1, radius1, sphereC1, kA, kD, kS);
	sphere sphere2(sphereV2, radius2, sphereC2, kA, kD, kS);
	sphere sphere3(sphereV3, radius3, sphereC3, kA, kD, kS);
	sphere sphere4(sphereV4, radius4, sphereC4, kA, kD, kS);
	sphere sphere5(sphereV5, radius5, sphereC5, kA, kD, kS);
	sphere sphere6(sphereV6, radius6, sphereC6, kA, kD, kS);
	sphere sphere7(sphereV7, radius7, sphereC7, kA, kD, kS);
	sphere sphere8(sphereV8, radius8, sphereC8, kA, kD, kS);
	sphere sphere9(sphereV9, radius9, sphereC9, kA, kD, kS);
	sphere sphere10(sphereV10, radius10, sphereC10, kA, kD, kS);
	sphere sphere11(sphereV11, radius11, sphereC11, kA, kD, kS);
	sphere sphere12(sphereV12, radius12, sphereC12, kA, kD, kS);
	sphere sphere13(sphereV13, radius13, sphereC13, kA, kD, kS);
	sphere sphere14(sphereV14, radius14, sphereC14, kA, kD, kS);
	sphere sphere15(sphereV15, radius15, sphereC15, kA, kD, kS);
	sphere sphere16(sphereV16, radius16, sphereC16, kA, kD, kS);
	sphere sphere17(sphereV17, radius17, sphereC17, kA, kD, kS);
	sphere sphere18(sphereV18, radius18, sphereC18, kA, kD, kS);
	sphere sphere19(sphereV19, radius19, sphereC19, kA, kD, kS);
	sphere sphere20(sphereV20, radius20, sphereC20, kA, kD, kS);
	spheres.push_back(sphere1);
	spheres.push_back(sphere2);
	spheres.push_back(sphere3);
	spheres.push_back(sphere4);
	spheres.push_back(sphere5);
	spheres.push_back(sphere6);
	spheres.push_back(sphere7);
	spheres.push_back(sphere8);
	spheres.push_back(sphere9);
	spheres.push_back(sphere10);
	spheres.push_back(sphere11);
	spheres.push_back(sphere12);
	spheres.push_back(sphere13);
	spheres.push_back(sphere14);
	spheres.push_back(sphere15);
	spheres.push_back(sphere16);
	spheres.push_back(sphere17);
	spheres.push_back(sphere18);
	spheres.push_back(sphere19);
	spheres.push_back(sphere20);

	//Ground plane and walls
	vec normal1(-5, 0, -4), point1(1300, 500, 500), groundColor(230, 182, 159);
	plane wall1(normal1, point1, groundColor, kA, kD, kS);
	vec normal2(5, 0, -4), point2(-1000, 500, 500);
	plane wall2(normal2, point2, groundColor, kA, kD, kS);
	vec normal3(0, 0, 1), point3(0, 0, 0);
	plane wall3(normal3, point3, groundColor, kA, kD, kS);
	vec normal(0, 5, -1), point(0, -60, 120);
	plane groundPlane(normal, point, groundColor, kA, kD, kS);
	std::vector<plane> planes{ groundPlane, wall1, wall2, wall3 };

	//Image properties
	const int antiAliasing = 3;
	const int H = 1080;  // Height of the image
	const int W = 1920;  // Width of the image
	std::ofstream image = createBmpFile(H, W);
	//Some settings for BMP file
	unsigned char bmpPad[3] = { 0, 0, 0 };
	const int paddingAmount = ((4 - (W * 3) % 4) % 4);

	//Screen, eye and light source properties
	double screenHeight = 100;          //Assuming same ratio for screen width too.
	vec upperLeftCorner(50, 50, 100);   //Position of the upper left corner of the screen
	double ratio = H / screenHeight;    //Ratio in pixels per coordinate unit, assuming delta x = delta y
	vec eye(0, 0, 0);                   //Position of the eye
	vec lightPosition(400, 400, 400);   //Position of the light source
	double I_a = 1, I_source = 585000;  //Ambient and light source intensities
	lightSource light(lightPosition, I_a, I_source);

	//Ray tracing
	for (int y = H; y >= 1; y--) {  //Rows of pixels, bottom to top (weird property of BMP)
		if (y % 50 == 0) {
			std::cerr << "\rLines remaining: " << y << " " << std::flush;
		}
		for (int x = 1; x <= W; x++) {  //Columns of pixels, left to right
			color pixelColor;
			//Anti-aliasing, random samples per pixel
			for (int j = 1; j <= antiAliasing; j++) {
				float dx = (float)rand() / RAND_MAX, dy = (float)rand() / RAND_MAX;
				vec tempPixel((-x + dx) / ratio, (-y + dy) / ratio, 0);                             //Position of the pixel relative to upper left corner.
				vec pixel = upperLeftCorner + tempPixel;                                            //Absolute position of the pixel.
				ray currentRay(eye, pixel - eye);                                                   //Create a ray from eye to center of the pixel.
				pixelColor = pixelColor + findColor(currentRay, spheres, planes, maxDepth, light);  //Check the intersection color
			}
			writeColor(pixelColor / antiAliasing, image);
		}
		image.write(reinterpret_cast<char*>(bmpPad), paddingAmount);  //Something about BMP
	}
	image.close();
	std::cin.get();
	return 0;
}