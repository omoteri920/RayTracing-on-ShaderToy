#define MAXX 10000000.0

float seed = 0.0;


vec3 bgColor(vec3 rayDir) {
    float u =  0.5*(1.0 + rayDir[1]);
    //return vec3(0.05, 0.05, 0.2);
    return u*vec3(0.7, 0.8, 0.9) + (1.0-u)*vec3(0.05, 0.05, 0.2);
}
    
// random number between 0 and 1
float random() {
    return fract(sin(seed++)*43758.5453123);
}

// a Light is defined by a location and a color
struct Light {
    vec3 location;
    vec3 color;
};

// Sphere is defined by a center and radius and material: color
struct Sphere {
	float radius;
	vec3 center;
    vec3 color;
};

// Ray is define by an origin point and a direction vector
struct Ray {
    vec3 origin;
    vec3 direction;
};
 
struct Intersection {
    int obj;
    float t;
};

const int NUM_SPHERE = 4;
Sphere spheres[NUM_SPHERE];
Light  lights[1];

// Intersection code for Ray-Sphere    
float raySphereIntersect(in Ray ray, in Sphere sphere) {
    
    vec3 rayToSphere = ray.origin - sphere.center;
    float b = dot(rayToSphere, ray.direction);
    float c = dot(rayToSphere, rayToSphere) - (sphere.radius * sphere.radius);
	float disc = b*b - c;
    float t;
    if (disc > 0.0) {
        t = -b - sqrt(disc);
        if (t > 0.00001) {
            return t;
        }
        t = -b + sqrt(disc);
        if (t > 0.00001) {
            return t;
        }  
    }
    return MAXX;
}

// Traverses the entire scene and 
// returns the objectID and the intersection point
Intersection intersectAllObjects(Ray ray) {
    float minT = MAXX;
    int iSphere = -1;
    
    for (int i=0; i < NUM_SPHERE; i++) {
       Sphere sphere = spheres[i];
       
       float t = raySphereIntersect(ray, sphere);
         
       if (t < minT && t >= 0.001) {
           // keep track of the closest sphere and intersection
           iSphere = i;
           minT = t;
       }
   }
   
   return Intersection(iSphere, minT);
}
  

// create 4 spheres at different locations in different colors
void makeScene(int f) {
    //spheres[0] = Sphere(1.0, vec3(0.2, 0, -7), vec3(1, 1, 0));
    //spheres[1] = Sphere(0.3, vec3(-0.4, -0.5, -2.5), vec3(1, 0, 0));
    spheres[2] = Sphere(0.4, vec3(0, 0, -2), vec3(0, 1, 0));
    //spheres[3] = Sphere(0.8, vec3(1., 1, -5.5), vec3(0, 0, 1));
    lights[0] = Light(vec3(0.0, 0.0, 0.0), vec3(1, 1, 1));
}

// samples per pixel
const int numSamples=4;



/* The main entry point:
   * This is called for every pixel on the screen 
*/
void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
   
   // fragCoord ranges from 
   //   in x: 0.5 to iResolution.x-0.5
   //   in y: 0.5 to iResolution.y-0.5
   // pixel (0,0) is at the bottom left corner
  
   makeScene(iFrame);
   
   //vec3 rayOrigin = vec3(0.0, 0.0, 0.0);
   
   float screenDepth = -2.0;
   
   float width = iResolution.x;
   float height = iResolution.y; 
   
   vec3 samp = vec3(0, 0, 0);
   seed = 0.0;
   for (int i=0; i<1*numSamples; i++) {
       float x = fragCoord.x + random() - 0.5;
       float y = fragCoord.y + random() - 0.5;
       
   
       // map (0.5, w-0.5) to (-1, 1)
       // and (0.5, h-0.5) to (-1, 1)
       x = (x/width)*2.0 - 1.0;
       y = (y/height)*2.0 - 1.0;
       
       // account for the non-square window
       float aspectRatio = width/height;
       y = y/aspectRatio;
       
       // ray Origin for is at (0, 0, 0)
       vec3 rayOrigin = vec3(0.0, 0.0, 0.0);
       // normalized ray direction
       vec3 rayDirection = normalize(vec3(x, y, screenDepth));
       
       Ray ray = Ray(rayOrigin, rayDirection);
       vec3 col = vec3(1.0, 1.0, 1.0);
       for (int k=0; k<100; k++) { // 100 bounces
       
       // traverse the scene (all spheres) and find the 
       // closest intersected object and intersection point
       Intersection intersection = intersectAllObjects(ray);
       
       int iSphere = intersection.obj;
       float minT = intersection.t;
       
       Sphere sphere;
       
       if (iSphere > -1) { 
           // to get around iSphere not being constant
           for (int i=0; i<NUM_SPHERE; i++) {
               if (i==iSphere) {
                   sphere = spheres[i];
                   break;
               }
           }
           
           // hit coordinates
           vec3 hit = ray.origin + minT*ray.direction;
           // normal at the point of ray-sphere intersection
           vec3 hitPointNormal = normalize(hit-sphere.center);
           vec3 xnorm = hitPointNormal;
           vec3 reflDir = ray.direction - 2.0*xnorm*dot(ray.direction, xnorm);
           ray.origin = hit;
           ray.direction = reflDir;
           /****
           // vector from eye to the intersection point
           vec3 hitPointToEye = normalize(vec3(0, 0, 0)-hit);
           // cosine of the angle between ray and nornal
           float angle = dot(hitPointNormal, hitPointToEye);
           // use the cosine of the angle to modulate color for 
           // a simple diffuse shading effect
           samp = samp + angle*sphere.color;
           ***/
           //samp = sphere.color;
           //col = col*sphere.color;        
       } 
       else {
           samp = samp + col*bgColor(ray.direction);
           break;
       }
   }
   }
   // average all the samples per pixel
   fragColor = vec4(samp/float(numSamples), 1.0);

}