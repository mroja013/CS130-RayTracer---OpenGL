#include "reflective_shader.h"
#include "ray.h"
#include "render_world.h"

vec3 Reflective_Shader::
Shade_Surface(const Ray& ray,const vec3& intersection_point,
    const vec3& normal,int recursion_depth) const
{
	vec3 final_color, shader_color, reflecteD;
	
	shader_color = shader->Shade_Surface(ray, intersection_point, normal, recursion_depth);

	reflecteD = (ray.direction - (2*dot(ray.direction,normal)*normal)).normalized();
	Ray reflected_ray = Ray(intersection_point, reflecteD);

	final_color = world.Cast_Ray(reflected_ray, recursion_depth+1);

	return final_color*reflectivity+(1-reflectivity)*shader_color;
}