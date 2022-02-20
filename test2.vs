#define MAXX 10000000.0

float seed = 0.0;

struct Material {
    vec3 baseColor;      
    float metallic;      
    float subsurface;    
    float specular;      
    float roughness;     
    float specularTint;  
    float anisotropic;   
    float sheen;         
    float sheenTint;     
    float clearcoat;     
    float clearcoatGloss;
};

const float PI = 3.14159265358979323846;

float sqr(float x) { return x*x; }

float SchlickFresnel(float u)
{
    float m = clamp(1.0-u, 0.0, 1.0);
    float m2 = m*m;
    return m2*m2*m;
}

float GTR1(float NdotH, float a)
{
    if (a >= 1.0) return 1.0/PI;
    float a2 = a*a;
    float t = 1.0 + (a2-1.0)*NdotH*NdotH;
    return (a2-1.0) / (PI*log(a2)*t);
}

float GTR2(float NdotH, float a)
{
    float a2 = a*a;
    float t = 1.0 + (a2-1.0)*NdotH*NdotH;
    return a2 / (PI * t*t);
}

float GTR2_aniso(float NdotH, float HdotX, float HdotY, float ax, float ay)
{
    return 1.0 / (PI * ax*ay * sqr( sqr(HdotX/ax) + sqr(HdotY/ay) + NdotH*NdotH ));
}

float smithG_GGX(float NdotV, float alphaG)
{
    float a = alphaG*alphaG;
    float b = NdotV*NdotV;
    return 1.0 / (NdotV + sqrt(a + b - a*b));
}

float smithG_GGX_aniso(float NdotV, float VdotX, float VdotY, float ax, float ay)
{
    return 1.0 / (NdotV + sqrt( sqr(VdotX*ax) + sqr(VdotY*ay) + sqr(NdotV) ));
}

vec3 mon2lin(vec3 x)
{
    return vec3(pow(x[0], 2.2), pow(x[1], 2.2), pow(x[2], 2.2));
}


vec3 BRDF( vec3 L, vec3 V, vec3 N, vec3 X, vec3 Y, Material m){
    float NdotL = clamp(L.z, 0.0, 1.0);
    float NdotV = clamp(V.z, 0.0, 1.0);
    if (NdotL < 0.0 || NdotV < 0.0) return vec3(0.0, 0.0, 0.0);

    vec3 H = normalize(L+V);
    float NdotH = clamp(H.z, 0.0, 1.0);
    float LdotH = clamp(dot(L, H), 0.0, 1.0);
    
    float HdotX = H.x;
    float HdotY = H.y; 
    

    vec3 Cdlin = mon2lin(m.baseColor);
    Cdlin = m.baseColor;
    float Cdlum = .3*Cdlin[0] + .6*Cdlin[1]  + .1*Cdlin[2]; 
    vec3 Ctint = Cdlum > 0.0 ? Cdlin/Cdlum : vec3(1);
    vec3 Cspec0 = mix(m.specular*.08*mix(vec3(1), Ctint, m.specularTint), Cdlin, m.metallic);
    vec3 Csheen = mix(vec3(1), Ctint, m.sheenTint);

    float FL = SchlickFresnel(NdotL), FV = SchlickFresnel(NdotV);
    float Fd90 = 0.5 + 2.0 * LdotH*LdotH * m.roughness;
    float Fd = mix(1.0, Fd90, FL) * mix(1.0, Fd90, FV);
    
    float fd_90_minus_1 = 2.0 * LdotH * LdotH * m.roughness - 0.5;
    
    Fd  = (1.0 + fd_90_minus_1 * pow(1.0 - NdotL, 5.0))
        * (1.0 + fd_90_minus_1 * pow(1.0 - NdotV, 5.0));

    float Fss90 = LdotH*LdotH*m.roughness;
    float Fss = mix(1.0, Fss90, FL) * mix(1.0, Fss90, FV);
    float ss = 1.25 * (Fss * (1.0 / (NdotL + NdotV) - .5) + .5);

    // specular
    float aspect = sqrt(1.0-m.anisotropic*.9);
    float ax = max(.001, sqr(m.roughness)/aspect);
    float ay = max(.001, sqr(m.roughness)*aspect);
    float Ds = GTR2_aniso(NdotH, HdotX, HdotY, ax, ay);
    float FH = SchlickFresnel(LdotH);
    vec3 Fs = mix(Cspec0, vec3(1), FH);
    float Gs;
    Gs  = smithG_GGX_aniso(NdotL, dot(L, X), dot(L, Y), ax, ay);
    Gs *= smithG_GGX_aniso(NdotV, dot(V, X), dot(V, Y), ax, ay);
    Gs = 1.0; 
    // sheen
    vec3 Fsheen = FH * m.sheen * Csheen;

    float Dr = GTR1(NdotH, mix(.1,.001,m.clearcoatGloss));
    float Fr = mix(.04, 1.0, FH);
    float Gr = smithG_GGX(NdotL, .25) * smithG_GGX(NdotV, .25);
    
    return ((1.0/PI) * mix(Fd, ss, m.subsurface)*Cdlin + Fsheen)
        * (1.0-m.metallic)
        + Gs*Fs*Ds + .25*m.clearcoat*Gr*Fr*Dr;
}


vec3 bgColor(vec3 rayDir) {
    float u =  0.5*(1.0 + rayDir[1]);
    //return vec3(0.05, 0.05, 0.2);
    return u*vec3(0.7, 0.8, 0.9) + (1.0-u)*vec3(0.05, 0.05, 0.2);
}
    
// random number between 0 and 1
float random() {
    return fract(sin(seed++)*43758.5453123);
}

struct Light {
    vec3 location;
    vec3 color;
    float radius;
};

struct Sphere {
	float radius;
	vec3 center;
    vec3 color;
};

struct Ray {
    vec3 origin;
    vec3 direction;
};
 
struct Intersection {
    int obj;
    float t;
};

const int NUM_SPHERE = 10;
Sphere spheres[NUM_SPHERE];
Material materials[NUM_SPHERE];
Light  lights[1];


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
    int count = 0;
    float radius = 0.25;
    float x = -0.5;
    float y = -0.2; 
    float z = -3.5;
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            spheres[count] = Sphere(radius,vec3(x, y, z), vec3(0.0));
            materials[count].baseColor = vec3(random(), random(), random());
            x = x + (radius+0.5);
            count++;
        }
        x = -0.5;
        z = z - (radius+0.5);
        y = y + 0.3;
    }
    //Ground
    spheres[9] = Sphere(50.0, vec3(0.0, -10.0, -70.0), vec3(1.0, 0.3, 0.3));

    //Light
    lights[0] = Light(vec3(.75, 1., 1.), vec3(1, 1, 1), 0.2);

    for(int i = 0; i < 9; i++){
        materials[i].baseColor = vec3(random()*.5, random()*.5, random()*.2);
        materials[i].metallic = 0.3; 
        materials[i].subsurface = 0.0;
        materials[i].specular = 0.2;
        materials[i].roughness = 0.4;
        materials[i].specularTint = 0.0;
        materials[i].anisotropic = 0.0;
        materials[i].sheen = 0.0;
        materials[i].sheenTint = 0.0;
        materials[i].clearcoat = 0.0;
        materials[i].clearcoatGloss = 0.0;
    }
    materials[9].baseColor = vec3(0.25, 0.25, 0.25);
}

// samples per pixel
const int numSamples=4;

void convertToTangentSpace(vec3 toLight, vec3 toView, vec3 hitPoint, out vec3 toLightTS, out vec3 toViewTS, out vec3 nTS) {
    vec3 t = normalize(dFdx(hitPoint));
    vec3 b = normalize(dFdy(hitPoint));
    vec3 n = normalize(cross(t, b));
    mat3 xformMatrix = transpose(mat3(t, b, n));
    
    toLightTS = xformMatrix * toLight;
    toViewTS = xformMatrix * toView;
    nTS = n;
}

float checkLightVisibility(in Light light, vec3 hitPoint, vec3 hitNormal) {
    float visible = 0.0;
    vec3 lightDir = normalize(light.location - hitPoint);
    Ray ray;
    ray.origin = hitPoint+lightDir*0.01;
    ray.direction = lightDir;
    Intersection intersection = intersectAllObjects(ray);
    int iSphere = intersection.obj;
    if (iSphere == -1) { 
        visible = 1.0;
    }
    return visible;
}


void mainImage(out vec4 fragColor, in vec2 fragCoord){
   makeScene(iFrame);
   
   vec3 rayOrigin = vec3(0.0, 0.0, 0.0);
   
   float screenDepth = -2.0;
   
   float width = iResolution.x;
   float height = iResolution.y; 
   
   vec3 samp = vec3(0, 0, 0);
   seed = 0.0;
   for (int i=0; i<1*numSamples; i++) {
       float x = fragCoord.x + random() - 0.5;
       float y = fragCoord.y + random() - 0.5;

       x = (x/width)*2.0 - 1.0;
       y = (y/height)*2.0 - 1.0;
       
       float aspectRatio = width/height;
       y = y/aspectRatio;
              
       vec3 rayDirection = normalize(vec3(x, y, screenDepth));
       
       Ray ray = Ray(rayOrigin, rayDirection);
              
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
           
           // hit coordinates
           vec3 hit = ray.origin + minT*ray.direction;
           // normal at the point of ray-sphere intersection
           vec3 hitPointNormal = normalize(hit-sphere.center);         
           // vector from intersection to light
           vec3 hitPointToLight = normalize(lights[0].location-hit);
           vec3 hitPointToView = ray.origin-hit;
           
           vec3 toViewTS;
           vec3 toLightTS;
           vec3 nTS;
           convertToTangentSpace(hitPointToLight, hitPointToView, hit, toLightTS, toViewTS, nTS); 
           vec3 h_ts 		= normalize(toLightTS + toViewTS);
           float dot_nl 	= clamp(toLightTS.z, 0.0, 1.0);
           float dot_nv 	= clamp(toViewTS.z, 0.0, 1.0);
           float dot_nh	= clamp(h_ts.z, 0.0, 1.0);
           float dot_lh 	= clamp(dot(toLightTS, h_ts), 0.0, 1.0);
           float dot_ht	= h_ts.x;
           float dot_hb	= h_ts.y;
           Material m = materials[iSphere];
           vec3 wi;
           float lightPdf;
  
           float lightVisible = checkLightVisibility(lights[0], hit, hitPointNormal);

           vec3 brdf = BRDF(toLightTS, toViewTS, nTS, vec3(0.0, 0.0, 0.0), vec3(0.0, 0.0, 0.0), m);
          
           vec3 col = brdf*dot_nl*lightVisible;
           vec3 toneMappedColor = col * (1.0 / (col + 1.0));
           float gamma = 1.0/2.2;
           vec3 finalColor = vec3(pow(toneMappedColor.x, gamma), 
                                  pow(toneMappedColor.y, gamma), 
                                  pow(toneMappedColor.z, gamma));
           samp = samp + finalColor; 
       } 
       else {
           samp = samp + bgColor(ray.direction);

       }
   }

   fragColor = vec4(samp/float(numSamples), 1.0);
}