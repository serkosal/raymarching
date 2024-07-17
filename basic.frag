#version 330 core

#define MIN_SURFACE 0.005
#define MAX_DIST 1000.0
#define MAX_STEPS 1000

#define MAX_OBJECT_COUNT 20
#define MAX_POINT_LIGHTS_COUNT 10

/***********************************************************************************
/**                                                                               **
/**    Structs block                                                              **
/**                                                                               **
/**********************************************************************************/

struct Sphere {
    vec3 color;
    vec3 pos;
    float radii;
};

struct Box {
    vec3 color;
    vec3 pos;
    vec3 size;
};

struct Torus {
    vec3 color;
    vec3 pos;
    vec2 size;
};

struct Inf_Cylinder {
    vec3 color;
    vec3 pos;
    vec3 size;
};

struct Cone {
    vec3 color;
    vec3 pos;
    vec2 base_center;
    float h;
};

struct Inf_Cone {
    vec3 color;
    vec3 pos;
    vec2 base_center;
};

struct Plane {
    vec3 color;

    vec3 normal;
    float h;
};

struct Hex_Prism {
    vec3 color;
    vec3 pos;
    vec2 h;
};

struct Tri_Prism {
    vec3 color;
    vec3 pos;
    vec2 h;
};

struct Capsule {
    vec3 color;
    vec3 pos;
    vec3 a;
    vec3 b;
    float r;
};

struct Capped_Cylinder {
    vec3 color;
    vec3 pos;
    vec2 h;
};

struct Point_Light {
    vec3 pos;

    float c, l, q;
    vec3 color;
};

/***********************************************************************************
/**                                                                               **
/**    Uniform, in/out variables block                                            **
/**                                                                               **
/**********************************************************************************/

uniform vec2 iResolution;
uniform float iTime;

uniform int spheres_num;
uniform Sphere spheres[MAX_OBJECT_COUNT];

uniform int boxes_num;
uniform Box boxes[MAX_OBJECT_COUNT];

uniform int planes_num;
uniform Plane planes[MAX_OBJECT_COUNT];

//uniform int toruses_num;
//uniform int inf_cylinders_num;
//uniform int cones_num;
//uniform int inf_cones_num;
//uniform int planes_num;
//uniform int hex_prisms_num;
//uniform int tri_prisms_num;
//uniform int capsules_num;
//uniform int capped_cylinders_num;

struct Camera
{
    vec3 pos;
    mat3 rot;
};

uniform Camera camera;

uniform int point_lights_num;
uniform Point_Light point_lights[MAX_POINT_LIGHTS_COUNT];

out vec4 frag_Color;

/***********************************************************************************
/**                                                                               **
/**    SDF declarations block                                                     **
/**                                                                               **
/**********************************************************************************/

float SDF_sphere(vec3 p, vec3 c, float r);
float SDF_box( vec3 p, vec3 size );
float SDF_torus( vec3 p, vec2 t );
float SDF_inf_cylinder( vec3 p, vec3 c );
float SDF_cone( vec3 p, vec2 c, float h );
float SDF_inf_cone( vec3 p, vec2 c );
float SDF_plane( vec3 p, vec3 n, float h );
float SDF_hex_prism( vec3 p, vec2 h );
float SDF_tri_prism( vec3 p, vec2 h );
float SDF_capsule( vec3 p, vec3 a, vec3 b, float r );
float SDF_caped_cylinder(vec3 p, float h, float r);

/***********************************************************************************
/**                                                                               **
/**    Rendering functions block                                                  **
/**                                                                               **
/**********************************************************************************/


float map_the_world(vec3 p, vec3 cam_pos)
{
    float dist = 1.0 / 0.0;
    //int object_type;
    //int object_num;

    float d;
    for (int i = 0; i < planes_num; i++)
    {
        d = SDF_plane(p, planes[i].normal, planes[i].h);
        if (d < dist)
        {
            dist = d;
        }
    }
    for (int i = 0; i < spheres_num; i++)
    {
        float d = SDF_sphere(p, spheres[i].pos, spheres[i].radii);
        if (d < dist)
        {
            dist = d;
        }
    }
    for (int i = 0; i < boxes_num; i++)
    {
        d = SDF_box(p - boxes[i].pos, boxes[i].size);
        if (d < dist)
        {
            dist = d;
        }
    }

    

    //mat3 rot = mat3(1.0); 
    //float phi = 0.25 * iTime;
    //rot[0] = vec3(cos(phi), 0.0, sin(phi));
    //rot[1] = vec3(0.0, 1.0, 0.0);
    //rot[2] = vec3(-sin(phi), 0.0, cos(phi));
    //p = rot * p;

    return dist;
}

vec3 calculate_normal(vec3 p, vec3 ro)
{
    const vec3 small_step = vec3(0.001, 0.0, 0.0);

    float gradient_x = map_the_world(p + small_step.xyy, ro) - map_the_world(p - small_step.xyy, ro);
    float gradient_y = map_the_world(p + small_step.yxy, ro) - map_the_world(p - small_step.yxy, ro);
    float gradient_z = map_the_world(p + small_step.yyx, ro) - map_the_world(p - small_step.yyx, ro);

    vec3 normal = vec3(gradient_x, gradient_y, gradient_z);

    return normalize(normal);
}

float shadow(vec3 lPos, vec3 lDir, float pDist)
{
    vec3 l = lPos;
    for (int i = 0; i < 50; i++)
    {
        float h = map_the_world(l, vec3(0.));
        if (h > MAX_DIST || h > pDist)
            break;
            //return 1.0;
        if (h < MIN_SURFACE)
        {
            if (abs(distance(l, lPos) - pDist) > 0.1)
                return 0.0;
            else return 1.0;
        }
        
        l += lDir * h;
    }

    return 1.0;
}

vec3 getLight(vec3 p, vec3 ro)
{
    vec3 n = calculate_normal(p, ro);

    float diffuse_max = 0.f;
    vec3 color;

    for (int i = 0; i < point_lights_num; ++i)
    {
        vec3 lPos = point_lights[i].pos;
        vec3 l = normalize(lPos - p);
        float pDist = distance(lPos, p);

        diffuse_max = max(max(dot(n, l), 0.01), diffuse_max);

        color += min(shadow(lPos, -l, pDist), diffuse_max) * point_lights[i].color 
        / (point_lights[i].c + (point_lights[i].l + point_lights[i].q * pDist) * pDist);
    }
    
    //return n * 0.5 + 0.5;
    
    return color;
}

vec3 rayMarch(vec3 ro, vec3 rd)
{
    vec3 r = ro;
    for (int i = 0; i < MAX_STEPS; i++)
    {
        float h = map_the_world(r, ro);
        if (h < MIN_SURFACE)
            return getLight(r, ro);
        
        if (h > MAX_DIST)
            break;
            
        r += rd*h;
    }
    
    return vec3(0.0);
}

/***********************************************************************************
/**                                                                               **
/**    Main function block                                                        **
/**                                                                               **
/**********************************************************************************/

void main()
{
    float zoom = 1.0;

	vec2 uv = gl_FragCoord.xy / iResolution; // uv in [0.f, 1.f]
	uv = uv * 2.0 - 1.0; // uv in [-1.f, 1.f]

    float aspectRatio = iResolution.x / iResolution.y;
    uv.x *= aspectRatio;


    //vec3 rd = camera.trans * normalize(vec3(uv.x, uv.y, zoom));
    vec3 rd = camera.rot * normalize(vec3(uv.x, uv.y, zoom));

    // Output to screen
    frag_Color = vec4(rayMarch(camera.pos, normalize(rd)), 1.0);
}

/***********************************************************************************
/**                                                                               **
/**    Signed-distance functions block                                            **
/**                                                                               **
/**********************************************************************************/

float SDF_sphere(vec3 p, vec3 c, float r)
{
    return distance(p, c.xyz) - r;
}

float SDF_box( vec3 p, vec3 size )
{
    vec3 d = abs(p) - size;
    return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
}

float SDF_torus( vec3 p, vec2 t )
{
    vec2 q = vec2(length(p.xz)-t.x,p.y);
    return length(q)-t.y;
}

float SDF_inf_cylinder( vec3 p, vec3 c )
{
  return length(p.xz-c.xy)-c.z;
}

float SDF_cone( vec3 p, vec2 c, float h )
{
  // c is the sin/cos of the angle, h is height
  // Alternatively pass q instead of (c,h),
  // which is the point at the base in 2D
  vec2 q = h*vec2(c.x/c.y,-1.0);
    
  vec2 w = vec2( length(p.xz), p.y );
  vec2 a = w - q*clamp( dot(w,q)/dot(q,q), 0.0, 1.0 );
  vec2 b = w - q*vec2( clamp( w.x/q.x, 0.0, 1.0 ), 1.0 );
  float k = sign( q.y );
  float d = min(dot( a, a ),dot(b, b));
  float s = max( k*(w.x*q.y-w.y*q.x),k*(w.y-q.y)  );
  return sqrt(d)*sign(s);
}

float SDF_inf_cone( vec3 p, vec2 c )
{
    // c is the sin/cos of the angle
    vec2 q = vec2( length(p.xz), -p.y );
    float d = length(q-c*max(dot(q,c), 0.0));
    return d * ((q.x*c.y-q.y*c.x<0.0)?-1.0:1.0);
}

float SDF_plane( vec3 p, vec3 n, float h )
{
  //n must be normalized
  return dot(p,n) + h;
  //return p.y;
}

float SDF_hex_prism( vec3 p, vec2 h )
{
  const vec3 k = vec3(-0.8660254, 0.5, 0.57735);
  p = abs(p);
  p.xy -= 2.0*min(dot(k.xy, p.xy), 0.0)*k.xy;
  vec2 d = vec2(
       length(p.xy-vec2(clamp(p.x,-k.z*h.x,k.z*h.x), h.x))*sign(p.y-h.x),
       p.z-h.y );
  return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

float SDF_tri_prism( vec3 p, vec2 h )
{
  vec3 q = abs(p);
  return max(q.z-h.y,max(q.x*0.866025+p.y*0.5,-p.y)-h.x*0.5);
}

float SDF_capsule( vec3 p, vec3 a, vec3 b, float r )
{
  vec3 pa = p - a, ba = b - a;
  float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
  return length( pa - ba*h ) - r;
}

float SDF_caped_cylinder(vec3 p, float h, float r)
{
    vec2 d = abs(vec2(length(p.xz),p.y)) - vec2(r,h);
    return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}
