#include "light.h"
#include "phong_shader.h"
#include "ray.h"
#include "render_world.h"
#include "object.h"
#include <vector>

vec3 Phong_Shader::
Shade_Surface(const Ray& ray,const vec3& intersection_point,
    const vec3& normal,int recursion_depth) const
{
    
    vec3 color, diffuse_component, specular_component, light_to_intersect, d, intensity;
    
    vec3 v, reflectVec, ambient_component;

	ambient_component = world.ambient_color*world.ambient_intensity*color_ambient;
	color = ambient_component;
    
    for(unsigned int i = 0; i < world.lights.size(); i++){

	//diffuse component calculations
        light_to_intersect = world.lights[i]->position - intersection_point;
        intensity = world.lights[i]->Emitted_Light(light_to_intersect);
        
        light_to_intersect = light_to_intersect.normalized();

	//forego remaining calculations if we have shadows enabled and this intersection is in shadow
	if(world.enable_shadows){

		Ray shine;
		shine.endpoint = world.lights[i]->position;
		shine.direction = (-light_to_intersect).normalized();
		Hit obscuring = world.Closest_Intersection(shine);
		double ogDist = (intersection_point - world.lights[i]->position).magnitude();
		double distance = (shine.Point(obscuring.dist) - world.lights[i]->position).magnitude();
		if (ogDist - distance > 0.00001) continue;

	}
        float temp = std::max(dot(light_to_intersect.normalized(), normal.normalized()), 0.0);
        
        diffuse_component = color_diffuse*intensity*temp;// add the diffuse values to color 

	//std::cout << "diffuse value: " << diffuse_component << std::endl;        
        
        //specular component calculations
        reflectVec = (- light_to_intersect + 2*(dot(light_to_intersect, normal))*normal).normalized();
        v = (world.camera.position - intersection_point).normalized();
        double highlight = (std::max(dot(reflectVec, v), 0.0));
        highlight = (pow(highlight, specular_power));
        specular_component = color_specular*intensity*highlight;
        
	//std::cout << "specular value: " << specular_component << std::endl;

        color += diffuse_component + specular_component;

//	std::cout << "total color value: " << color << std::endl;
        
        
        
    }
    
   // std::cout << "final return color: " << color << std::endl;    
    return color;
    
}
