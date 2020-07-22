#include "sphere.h"
#include "ray.h"

// Determine if the ray intersects with the sphere
Hit Sphere::Intersection(const Ray& ray, int part) const
{
    
    vec3 p = ray.endpoint - center;
    float alpha = pow(dot(ray.direction, p), 2) - (dot(ray.direction, ray.direction)) * (dot(p, p) - pow(radius, 2));
    if (alpha > 0){
    
      float t1 = -(dot(ray.direction, p)) + sqrt(alpha);
      float t2 = -(dot(ray.direction, p)) - sqrt(alpha);
      
      Hit intersect1, intersect2;
      
      if(t1 >= 0){
      
        intersect1.object = this;
        intersect1.dist = t1;
        
      }
      
      if (t2 >= 0){
      
        intersect2.object = this;
        intersect2.dist = t2;
      
      }
      
      if (intersect1.dist > intersect2.dist)
        return intersect2;
      
      return intersect1; 
      
    }
    Hit null;
    return null;
}

vec3 Sphere::Normal(const vec3& point, int part) const
{
    vec3 normal;
    normal = point - center;
    normal = normal.normalized();
    return normal;
}

Box Sphere::Bounding_Box(int part) const
{
    Box box;
    TODO; // calculate bounding box
    return box;
}
