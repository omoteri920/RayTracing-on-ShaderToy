#define MAXX 10000000.0

float seed = 0.0;
const float PI = 3.141592;


vec3 bgColor(vec3 rayDir) {
    float u =  0.5*(1.0 + rayDir[1]);
    return u*vec3(0.7, 0.8, 0.9) + (1.0-u)*vec3(0.05, 0.05, 0.2);
}
    
float random() {
    return fract(sin(seed++)*43758.5453123);
}

vec3 random_unit_vector(){
    float a = 2.*PI*random();
    float z = random()*2.-1.;
    float r = sqrt(1. - z*z);
    return vec3(r*cos(a), r*sin(a), z);
}

//material : 0→ランバート 1→スペキュラー 2→ライト(簡易)

struct Light {
    vec3 location;
    vec3 color;
};

struct Sphere_Light {
    float radius;
	vec3 center;
    vec3 color;
    int material;
};

struct Sphere {
	float radius;
	vec3 center;
    vec3 color;
    int material;
};

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
Sphere_Light  lights[1];

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


Intersection intersectAllObjects(Ray ray) {
    float minT = MAXX;
    int iSphere = -1;
    
    for (int i=0; i < NUM_SPHERE; i++) {
       Sphere sphere = spheres[i];
       
       float t = raySphereIntersect(ray, sphere);
         
       if (t < minT && t >= 0.001) {
           iSphere = i;
           minT = t;
       }
   }
   
   return Intersection(iSphere, minT);
} 

void makeScene(int f) {
    //spheres[2] = Sphere(1.0, vec3(0.2, 0, -7), vec3(1, 1, 0));
    spheres[0] = Sphere(0.3, vec3(.0, -.5, -2.), vec3(0, 1, 0), 0);
    spheres[1] = Sphere(0.4, vec3(.0, .6, -2.5), vec3(1, 0, 0), 1);
    spheres[2] = Sphere(0.2, vec3(-.5, .0, -2.2), vec3(1, 0, 0), 1);
    //spheres[3] = Sphere(0.8, vec3(1., 1, -5.5), vec3(0, 0, 1));
    //lights[0] = Light(vec3(0.0, 0.0, 0.0), vec3(1, 1, 1));
    spheres[3] = Sphere(0.5, vec3(1.5, .0, -2.), vec3(1.0, 1.0, 1.0), 2);
    //spheres[3] = Sphere(0.5, vec3(.0, 1.5, -2.), vec3(1.0, 1.0, 1.0), 2);
}

const int numSamples=50;


void mainImage(out vec4 fragColor, in vec2 fragCoord) {
  
   makeScene(iFrame);
   
   float screenDepth = -2.0;
   
   float width = iResolution.x;
   float height = iResolution.y; 
   
   vec3 samp = vec3(0, 0, 0);
   seed = 0.0;
   vec3 attenuation;
   //float attenuation_value = .5;
   for (int i=0; i<1*numSamples; i++) {
       float x = fragCoord.x + random() - 0.5;
       float y = fragCoord.y + random() - 0.5;
       
       x = (x/width)*2.0 - 1.0;
       y = (y/height)*2.0 - 1.0;
       
       float aspectRatio = width/height;
       y = y/aspectRatio;
       
       vec3 rayOrigin = vec3(0.0, 0.0, 0.0);
       
       vec3 rayDirection = normalize(vec3(x, y, screenDepth));
       
       Ray ray = Ray(rayOrigin, rayDirection);
       vec3 col = vec3(1.0, 1.0, 1.0);

       attenuation = vec3(.9, .9, .9);

       for (int k=0; k<100; k++) {
            Intersection intersection = intersectAllObjects(ray);
            
            int iSphere = intersection.obj;
            float minT = intersection.t;
            
            Sphere sphere;
            
            if (iSphere > -1) { 
                for (int i=0; i<NUM_SPHERE; i++) {
                    if (i==iSphere) {
                        sphere = spheres[i];
                        break;
                    }
                }

                if(sphere.material == 2){
                    samp = samp + attenuation * sphere.color;
                    break;
                }
            
                vec3 hit = ray.origin + minT*ray.direction;
                
                vec3 hitPointNormal = normalize(hit-sphere.center);
                vec3 xnorm = hitPointNormal;

                vec3 reflDir;
               
                if(sphere.material == 0) reflDir = hit + xnorm + random_unit_vector();
                else if(sphere.material == 1) reflDir = ray.direction - 2.0*xnorm*dot(ray.direction, xnorm) + 0.2 * random_unit_vector();
                else reflDir = hit + xnorm + random_unit_vector();

                ray.origin = hit;
                ray.direction = reflDir;
                attenuation = attenuation * sphere.color;
            } 
            else {
                samp = samp + attenuation * vec3(.1, .1, .1) * bgColor(ray.direction);
                break;
            }
        }
   }
   vec3 Color = (samp/float(numSamples));
   Color.x = sqrt(Color.x); Color.y = sqrt(Color.y); Color.z = sqrt(Color.z);
   fragColor = vec4(Color, 1.0);
}