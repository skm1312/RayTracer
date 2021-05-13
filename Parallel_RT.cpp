// Copyright (C) 2012  www.scratchapixel.com
// Public GNU license

// Modified by Krishna Murthy Surya Narayanan, 2020. 
// Last changed: 2020-07-31
/*
	 * Project:  Simple Ray Tracer
	 * License:  Public GNU
	 * Author:   www.scratchapixel.com
	 * Description: This is the work that i follwed for the implementation of my Serial code from which the following parallel code has been devised.
	 */
	 
/* A reference from the following work has been considered for this purpose as well

	/* File:     TinyKaboom.cpp
	 * Project:  TinyKaboom
	 * License:  Public GNU
	 * Author:   Dmitry V. Sokolov
	 * Date:     Jan-30-2019
	 * Description: The methods used are the sphere genration and rendering help.
	 */

#define _USE_MATH_DEFINES
#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "geometry.h"

struct Light {
    Light(const Vec3f &p, const float i) : position(p), intensity(i) {}
    Vec3f position;
    float intensity;
};

struct Surface_C {
    Surface_C(const float r, const Vec4f &a, const Vec3f &color, const float spec) : refractive_index(r), albedo(a), diffuse_color(color), specular_exponent(spec) {}
    Surface_C() : refractive_index(1), albedo(1,0,0,0), diffuse_color(), specular_exponent() {}
    float refractive_index;
    Vec4f albedo;
    Vec3f diffuse_color;
    float specular_exponent;
};

struct Sphere {
    Vec3f center;							//Center of the spehere
    float radius;							//Radius of the spehere
    Surface_C surface_C;					//Surface Color is explored here

    Sphere(const Vec3f &c, const float r, const Surface_C &m) : center(c), radius(r), surface_C(m) {}

    bool ray_intersect(const Vec3f &orig, const Vec3f &dir, float &t0) const {					// We identify the Ray-spehere intersection using euqations
        Vec3f L = center - orig;
        float tca = L*dir;
        float d2 = L*L - tca*tca;
        if (d2 > radius*radius) return false;
        float thc = sqrtf(radius*radius - d2);
        t0       = tca - thc;
        float t1 = tca + thc;
        if (t0 < 0) t0 = t1;
        if (t0 < 0) return false;
        return true;
    }
};

Vec3f reflect(const Vec3f &I, const Vec3f &N) {										// We recursively call for the reflect function, to obtain a more detailed view.
    return I - N*2.f*(I*N);
}


bool scene_intersect(const Vec3f &orig, const Vec3f &dir, const std::vector<Sphere> &spheres, Vec3f &hit, Vec3f &N, Surface_C &surface_C) {
    float spheres_dist = std::numeric_limits<float>::max();
    for (size_t i=0; i < spheres.size(); i++) {
        float dist_i;
        if (spheres[i].ray_intersect(orig, dir, dist_i) && dist_i < spheres_dist) {
            spheres_dist = dist_i;
            hit = orig + dir*dist_i;
            N = (hit - spheres[i].center).normalize();
            surface_C = spheres[i].surface_C;
        }
    }
	float plane_dist = std::numeric_limits<float>::max();
    return std::min(spheres_dist, plane_dist)<1000;
}

Vec3f cast_ray(const Vec3f &orig, const Vec3f &dir, const std::vector<Sphere> &spheres, const std::vector<Light> &lights, size_t depth=0) {		// The ray casting is done here and we deterimne if any scene has been intersected
    Vec3f point, N;
    Surface_C surface_C;

    if (depth>4 || !scene_intersect(orig, dir, spheres, point, N, surface_C)) {
        return Vec3f(0, 0, 0); // background color
    }

    Vec3f reflect_dir = reflect(dir, N).normalize();
    Vec3f reflect_orig = reflect_dir*N < 0 ? point - N*1e-3 : point + N*1e-3; 
    Vec3f reflect_color = cast_ray(reflect_orig, reflect_dir, spheres, lights, depth + 1);

    float diffuse_light_intensity = 0, specular_light_intensity = 0;
    for (size_t i=0; i<lights.size(); i++) {
        Vec3f light_dir      = (lights[i].position - point).normalize();
        float light_distance = (lights[i].position - point).norm();

        Vec3f shadow_orig = light_dir*N < 0 ? point - N*1e-3 : point + N*1e-3; 						// Verifying if shadow of the light is intersecting with any point
        Vec3f shadow_pt, shadow_N;
        Surface_C tmpmaterial;
        if (scene_intersect(shadow_orig, light_dir, spheres, shadow_pt, shadow_N, tmpmaterial) && (shadow_pt-shadow_orig).norm() < light_distance)
            continue;

        diffuse_light_intensity  += lights[i].intensity * std::max(0.f, light_dir*N);
        specular_light_intensity += powf(std::max(0.f, -reflect(-light_dir, N)*dir), surface_C.specular_exponent)*lights[i].intensity;
    }
    return surface_C.diffuse_color * diffuse_light_intensity * surface_C.albedo[0] + Vec3f(1., 1., 1.)*specular_light_intensity * surface_C.albedo[1] + reflect_color*surface_C.albedo[2];
}

	 
void render(const std::vector<Sphere> &spheres, const std::vector<Light> &lights) { 
	clock_t Tf;							
	Tf=clock();
	const int   width    = 1024;
    const int   height   = 768;
    const float fov      = M_PI/3.;
	int thrd			 = 7;					// Change the number of threads here.
	float Ts, Tcomp		 = 0;
	int W, count, Ds	 = 0;					// Variables for Computing the parallel run time
    std::vector<Vec3f> framebuffer(width*height);

    			
    for (size_t j = 0; j<height; j++) { 		// No of threads for parallel execution
        for (size_t i = 0; i<width; i++) {
            float dir_x =  (i + 0.5) -  width/2.;
            float dir_y = -(j + 0.5) + height/2.;    	
            float dir_z = -height/(2.*tan(fov/2.));
			count++;
            framebuffer[i+j*width] = cast_ray(Vec3f(0,0,0), Vec3f(dir_x, dir_y, dir_z).normalize(), spheres, lights);
        }
    }

    std::ofstream ofs; 									// The buffer keeps track of the file queue information
    ofs.open("./parallel.ppm",std::ios::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (size_t i = 0; i < height*width; ++i) {
        Vec3f &c = framebuffer[i];
        float max = std::max(c[0], std::max(c[1], c[2]));
        if (max>1) c = c*(1./max);
        for (size_t j = 0; j<3; j++) {
            ofs << (char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }
    ofs.close();
	
	Tf = clock()-Tf;								// Gives us the run time for rendering
	W = Ds + count;
	Ts = W * (((float)Tf/CLOCKS_PER_SEC)/100.00);	// Serial time in milli seconds
	Tcomp = Ts/ thrd;								// Parallel time in milli seconds
	
	std::cout<<"Parallel time complexity for "<<" "<<thrd<<" threads is"<< Tcomp<<" milliseconds";
}

int main() {
    Surface_C      red(1.0, Vec4f(0.6,  0.3, 0.1, 0.0), Vec3f(1.00, 0.32, 0.36),   50.);
    Surface_C 	  cyan(1.0, Vec4f(0.9,  0.1, 0.0, 0.0), Vec3f(0.65, 0.77, 0.97),   10.);
	Surface_C      violet(1.5, Vec4f(0.7,  0.5, 0.1, 0.8), Vec3f(0.42, 0, 1),  4.);

    std::vector<Sphere> spheres;
	spheres.push_back(Sphere(Vec3f( -3.0, 0, -20), 3, red));
	spheres.push_back(Sphere(Vec3f( 5.0,0, -25), 4, cyan));
	spheres.push_back(Sphere(Vec3f( 2.5,-4, -12), 2, violet));

    std::vector<Light>  lights;									// The settings mentioned in the above stated reference.
    lights.push_back(Light(Vec3f(-20, 20,  20), 1.5));
    lights.push_back(Light(Vec3f( 30, 50, -25), 1.8));
    lights.push_back(Light(Vec3f( 30, 20,  30), 1.7));

    render(spheres, lights);

    return 0;
}