// Author:
// Title:
// https://raytracing.github.io/
// https://github.com/RayTracing/raytracing.github.io/
// https://images.math.cnrs.fr/Un-arbre-pythagoricien.html
// https://www.josleys.com/
// https://cal.cs.umbc.edu/Courses/CMSC435-F15/Slides/raytrace.pdf
// https://www.cs.cornell.edu/courses/cs4620/2014fa/lectures/05rt-shading.pdf

// http://www.joelsornette.fr/Archives/exotypes/exotype26.pdf

// Ray Marching :
// https://iquilezles.org/articles/distfunctions/ 

// #ifdef GL_ES
// precision mediump float;
precision highp float;
// #endif

uniform sampler2D hokusai;
uniform sampler2D sky;
uniform sampler2D earth;
uniform sampler2D ground;
uniform sampler2D sun;
uniform sampler2D galaxy;
uniform sampler2D moon;
uniform sampler2D sagittarius_A;

uniform bool bases_vectors ;
uniform int scene ;


uniform vec2 u_resolution;
uniform vec2 u_mouse;
uniform float u_time;
uniform vec2 position;
uniform float light_height;
uniform int light_type; // 0: no light, 1: light+ambiente, 2: only light
uniform int draw_type ;

////////////////////////////////////////////////////////////////////////
//                                                                    //
//                            CONSTANT.glsl                           //
//                                                                    //
//                         Defining constants                         //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#define PI 3.1415926535897932384626433832795
#define TWO_PI 6.283185307179586476925286766559

#define COLOR_HOKUSAI vec3(0.0, 0.0, 0.01)
#define COLOR_EARTH vec3(0.0, 0.0, 0.02)
#define COLOR_SKY vec3(0.0, 0.0, 0.03)
#define COLOR_GROUND vec3(0, 0, 0.04)
#define COLOR_SUN vec3(0, 0, 0.05)
#define COLOR_GALAXY vec3(0, 0, 0.06)
#define COLOR_MOON vec3(0, 0, 0.07)
#define COLOR_SAGITTARIUS vec3(0, 0, 0.08)

#define COLOR_LATTICE_1 vec3(0,0.01,0)
#define COLOR_LATTICE_2 vec3(0,0.02,0)
#define COLOR_LATTICE_3 vec3(0,0.03,0)
#define COLOR_LATTICE_4 vec3(0,0.04,0)
#define COLOR_LATTICE_5 vec3(0,0.05,0)


#define TYPE_ERROR -1
#define TYPE_GROUPS 0
#define TYPE_SPHERE 1
#define TYPE_CYLINDER 2
#define TYPE_PLANE 3
#define TYPE_RECTANGLE 4
#define TYPE_CIRCLE 5
#define TYPE_TRIANGLE 6
#define TYPE_MANDELBULB 7

#define STRUCT_NULL -1
#define STRUCT_LIST 0
#define STRUCT_CUBE 1
#define STRUCT_CYLINDER 2
#define STRUCT_UNION 3
#define STRUCT_INTERSECTION 4
#define STRUCT_SMOOTH_U 5
#define STRUCT_CROPPED 6

// Structure definitions
struct Line {
    vec3 origin ;
    vec3 v ;
};

struct Sphere {
    vec3 center ;
    float radius ;
};

struct Cylinder {
    vec3 origin ;
    vec3 v ;
    float radius ;
};

struct Plane {
    vec3 origin ;
    vec3 u1 ;
    vec3 u2 ;
};


struct Object {
    Sphere sphere ;
    Plane plane ;
    Cylinder cylinder ;
    int type ;
    vec3 color ;
    float mirror ;
    float refrac ;
    bool is_light ;
};
#define GROUP_SIZE 10
struct Group {
    Object objects[GROUP_SIZE] ;
    int nb_objects ;
    int type ;
};

////////////////////////////////////////////////////////////////////////
//                                                                    //
//                       ObjectsInitializer.glsl                      //
//                                                                    //
//                 Functions for initialize an object                 //
//                                                                    //
////////////////////////////////////////////////////////////////////////

// Initialize an void object
Object null_Object(){
    return Object(
        Sphere(vec3(0.),0.),                      // sphere
        Plane(vec3(0.),vec3(0.),vec3(0.)),    // plane
        Cylinder(vec3(0.),vec3(0.),0.),         // cylinder
        TYPE_ERROR,                                 // type
        vec3(0.0),                                // color
        0.0,                                        // mirror
        1.0,                                        // refrac
        false                                       // is_light
    ) ;
}

// Initialize a sphere
Object new_Object(Sphere S, vec3 color, float mirror){
    Object O = null_Object();
    O.sphere = S;
    O.type = TYPE_SPHERE;
    O.color = color;
    O.mirror = mirror;
    return O;
}

// Initialize a cylinder
Object new_Object(Cylinder C, vec3 color, float mirror){
    Object O = null_Object();
    O.cylinder = C;
    O.type = TYPE_CYLINDER;
    O.color = color;
    O.mirror = mirror;
    return O;
}

// Initialize a plane
Object new_Object(Plane P, vec3 color, float mirror){
    Object O = null_Object();
    O.plane = P;
    O.type = TYPE_PLANE;
    O.color = color;
    O.mirror = mirror;
    return O;
}

// Initialize a rectangle
Object new_Rectangle(Plane P, vec3 color, float mirror){
    Object O = new_Object(P, color, mirror);
    O.type = TYPE_RECTANGLE;
    return O;
}
// Initialize a plane mirror
Object new_Plane(Plane P, vec3 color, float mirror){
    Object O = new_Object(P, color, mirror);
    O.type = TYPE_PLANE;
    return O;
}
// Initialize a circle / ellipse 
Object new_Circle(Plane P, vec3 color, float mirror){
    Object O = new_Object(P, color, mirror);
    O.type = TYPE_CIRCLE;
    return O;
}
// Initialize a plan and its type of rendering
Object new_typed_Plan(Plane P, vec3 color, float mirror, int type){
    Object O = new_Object(P, color, mirror);
    O.type = type;
    return O;
}
// Initialize a mandelbulb 
Object new_Bulb(vec3 pos, vec3 color, float mirror){
    Object O = null_Object();
    O.sphere = Sphere(pos,0);
    O.type = TYPE_MANDELBULB;
    O.color = color;
    O.mirror = mirror;
    return O;
}

// set the refraction of a object
Object Refractor(Object O, float refrac){
    O.refrac = refrac;
    return O;
}
// add light on an object
Object Light(Object O){
    O.is_light = true;
    return O ;
}


Group null_Group(){
    return Group(
        Object[10](
            null_Object(),null_Object(),null_Object(),null_Object(),null_Object(),
            null_Object(),null_Object(),null_Object(),null_Object(),null_Object()
        ),
        0,
        STRUCT_NULL
    );
}

Group Cube(vec3 pos,vec3 u1,vec3 u2,vec3 u3, vec3 col, float mirror){
    Group G = null_Group();
    G.type = STRUCT_CUBE;
    G.nb_objects = 6;
    G.objects[0] = new_Rectangle(Plane(pos,u1,u2),col,mirror);
    G.objects[1] = new_Rectangle(Plane(pos,u2,u3),col,mirror);
    G.objects[2] = new_Rectangle(Plane(pos,u3,u1),col,mirror);
    G.objects[3] = new_Rectangle(Plane(pos+u1+u2+u3,-u1,-u2),col,mirror);
    G.objects[4] = new_Rectangle(Plane(pos+u1+u2+u3,-u2,-u3),col,mirror);
    G.objects[5] = new_Rectangle(Plane(pos+u1+u2+u3,-u3,-u1),col,mirror);
    return G;
}

float norm2(vec3 v);
vec3 normalize(vec3 v);

Group Full_Cylinder(Cylinder C, vec3 col, float mirror){
    Group G = null_Group();
    G.type = STRUCT_CYLINDER;
    G.nb_objects = 3;

    vec3 v1 = cross(vec3(1,0,0),C.v);
    vec3 v2 = cross(vec3(0,1,0),C.v);
    if (norm2(cross(v1,v2))==0.0) v2 = cross(vec3(0,0,1),C.v);
    
    G.objects[0] = new_Object(C,col,mirror);
    G.objects[1] = new_Circle(
        Plane(C.origin,C.radius*normalize(v1),C.radius*normalize(v2)),
        col,mirror
    );
    G.objects[2] = new_Circle(
        Plane(C.origin+C.v,C.radius*normalize(v1),C.radius*normalize(v2)),
        col,mirror
    );
    return G;
}


////////////////////////////////////////////////////////////////////////
//                                                                    //
//                             Scenes.glsl                            //
//                                                                    //
//                  Creating the differents scenes                    //
//                     and placing objects in it                      //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#define nb_object 15

// Idea :
// List of objects
// List of links
// obj : [0,1,2,3,4,5]
// links : [union 0 1, inter 2 3 , ...]
// Calc : [o0, o1, o2, o3, o4, o5, union 0 1, inter 7 2, ...]
// Calc : Array[3*n]

// The array of all objects to render
Object objects[nb_object] = Object[nb_object](
    null_Object(),
    null_Object(),
    null_Object(),
    null_Object(),
    null_Object(),
    null_Object(),
    null_Object(),
    null_Object(),
    null_Object(),
    null_Object(),
    null_Object(),
    null_Object(),
    null_Object(),
    null_Object(),
    null_Object()
) ;

int index_obj = 0 ; 

// Push an object into objects, view as a chained list
void push_object(Object obj){
    // inout int index_obj, 
    objects[index_obj] = obj ;
    index_obj++ ;
}

// create and push the six faces of a cube
void push_cube(vec3 pos,vec3 u1,vec3 u2,vec3 u3, vec3 col, float mirror){
    push_object(new_Rectangle(Plane(pos,u1,u2),col,mirror));
    push_object(new_Rectangle(Plane(pos,u2,u3),col,mirror));
    push_object(new_Rectangle(Plane(pos,u3,u1),col,mirror));
    push_object(new_Rectangle(Plane(pos+u1+u2+u3,-u1,-u2),col,mirror));
    push_object(new_Rectangle(Plane(pos+u1+u2+u3,-u2,-u3),col,mirror));
    push_object(new_Rectangle(Plane(pos+u1+u2+u3,-u3,-u1),col,mirror));
}

// create and push a cylinder a his two extremity
void push_full_cylinder(Cylinder C, vec3 col, float mirror) {
    vec3 v1 = cross(vec3(1,0,0),C.v);
    vec3 v2 = cross(vec3(0,1,0),C.v);
    vec3 u = cross(v1,v2);
    if (dot(u,u)==0) v2 = cross(vec3(0,0,1),C.v);
    
    push_object(new_Object(C,col,mirror));
    push_object(new_Circle(
        Plane(C.origin,C.radius*normalize(v1),C.radius*normalize(v2)),
        col,mirror
    ));
    push_object(new_Circle(
        Plane(C.origin+C.v,C.radius*normalize(v1),C.radius*normalize(v2)),
        col,mirror
    ));

}

// Initialize the differente scene
void init_objects(){

    // Bases vectors
    if (bases_vectors){
        push_object(new_Object(Cylinder(vec3(0,0,0),vec3(1,0,0),.03),vec3(1,0,0),0.0));
        push_object(new_Object(Cylinder(vec3(0,0,0),vec3(0,1,0),.03),vec3(0,1,0),0.0));
        push_object(new_Object(Cylinder(vec3(0,0,0),vec3(0,0,1),.03),vec3(0,0,1),0.0));
    }


    // Scene 1 : somes objects
    if (scene == 1){
        push_object(Light(new_Object(
            //Sphere(vec3(0,light_height/2.5,0),1.2),
            Sphere(vec3(0,light_height/2.5,0),.8),
            COLOR_SUN,
            // vec3(1.0, 0.8157, 0.0),
            0.0
        )));
        push_object(new_Rectangle(
            Plane(vec3(-5,1,5),vec3(10,0,0),vec3(0,6,0)),
            COLOR_GALAXY, 0.0
        ));
        // push_object(new_Object(Sphere(vec3(-2,1.5,0),1),vec3(0.0, 1.0, 0.702),0.0));
        // push_object(new_Object(Sphere(vec3(2.5,3.5,-1),.4),vec3(0.6588, 0.0706, 0.6392),0.0));

        push_object(new_Object(Sphere(vec3(-2,1.5,0),1),COLOR_EARTH,0.0));
        push_object(new_Object(Sphere(vec3(2.5,3.5,-1),.4),COLOR_MOON,0.0));
        
        // push_object(new_Object(Cylinder(vec3(-4,1,3),vec3(2,3,-1),2),COLOR_SAGITTARIUS,0));
        // push_object(new_Object(Cylinder(vec3(0,0,0),vec3(0,1.01,0),1),COLOR_SKY,0));
    }
    

    // Scene 2 : Mirrors example
    if (scene == 2){
        push_object(Light(new_Object(Sphere(vec3(0,light_height/2+4,0),.2),COLOR_SUN,0.0)));
        push_object(new_Rectangle(
           Plane(vec3(-7.5,0,7),vec3(15,0,0),vec3(0,10,0)),
           //Plane(vec3(-7.5,0,7),vec3(15,0,0),vec3(0,10,-.5)),
           vec3(0.0, 0.0, 0.0),.5
        ));
        push_object(new_Object(Sphere(vec3(-2,1.5,0),1),COLOR_EARTH,0.0));
        push_object(new_Object(Sphere(vec3(2,2.2,1),1.8),vec3(1.0, 0.6, 0.0),.5));
        // push_object(new_Object(Cylinder(vec3(-4,1,1),vec3(2,3,-1),.5),vec3(0.0, 0.7, 1.0),.5));
        push_full_cylinder(
            Cylinder(vec3(-4,1,1),vec3(2,3,-1),.5),
            vec3(0.0, 0.7, 1.0),
            .5
        )
    }
    
    // Scene 3 : Refraction example
    if (scene == 3){

        // float refrac = light_height / 10.0 - 1. + 0.01 ;
        float refrac = exp(light_height/5) + 0.001;

        push_object(Refractor(new_Rectangle(
           Plane(vec3(0,1.5,5),vec3(10,0,0),vec3(0,5,0)),
           vec3(0,0,0),.9
        ), refrac));
        push_object(Refractor(new_Rectangle(
           Plane(vec3(0,1.5,5),vec3(0,0,10),vec3(0,5,0)),
           vec3(0,0,0),.9
        ), refrac));
        push_object(new_Object(Cylinder(vec3(-1,1,8),vec3(10,10,-2),.5),vec3(0.0, 0.7, 1.0),0));
    }

    // Scene 4 : A mandelbulb
    if (scene == 4){
        push_object(new_Bulb(vec3(-1,1,0),vec3(1),0)) ;
        vec3 col = vec3(0.3333, 0.7176, 0.5294);
        push_object(new_Object(Sphere(vec3(2,light_height,0),1.5),col,0.0));
    }

    // Scene 5 : A cube and intersection for ray marching
    if (scene == 5){
        vec3 col = vec3(0.3333, 0.7176, 0.5294);
        push_cube(vec3(-2,1,0),vec3(2,0,0),vec3(0,2,0),vec3(0,0,2),col,0.2);
        push_object(new_Rectangle(
           Plane(vec3(-7.5,0,7),vec3(15,0,0),vec3(0,10,0)),
           vec3(0.9569, 0.6627, 0.251),.5
        ));
        push_full_cylinder(
            Cylinder(vec3(4,1,1),vec3(-2,0,10),.5),
            vec3(0.0, 0.7, 1.0),.2
        );
    }

    // Scene 6 : A room with mirrors walls
    if (scene == 6){
        
        float z = 10 ;
        float dy = 10 ;
        float dx = 15 ;
        float dz = -10 ;

        push_object(new_Rectangle( // ahead
           Plane(vec3(-dx/2,0,z),vec3(dx,0,0),vec3(0,dy,0)),
           vec3(1.0, 0.0, 0.0),1
        ));
        push_object(new_Rectangle( // 
           Plane(vec3(-dx/2,0,z+dz),vec3(dx,0,0),vec3(0,dy,0)),
           vec3(1.0, 0.0, 0.0),1
        ));

        push_object(new_Rectangle( // top
           Plane(vec3(-dx/2,dy,z),vec3(dx,0,0),vec3(0,0,dz)),
           vec3(0.5, 0.5, 0.5),.0
        ));

        push_object(new_Rectangle( // left
           Plane(vec3(-dx/2,0,z),vec3(0,0,dz),vec3(0,dy,0)),
           vec3(1.0, 1.0, 0.0),0
        ));
        push_object(new_Rectangle( // right
           Plane(vec3(dx/2,0,z),vec3(0,0,dz),vec3(0,dy,0)),
           vec3(0.0, 0.0, 1.0),0
        ));

        push_object(new_Object(Sphere(vec3(5,2,z-4),1.5),vec3(0,1,0),0.6));
        push_object(new_Object(Cylinder(vec3(-4,1,z-5),vec3(2,7,-2),.5),vec3(0.0, 0.7, 1.0),0.5));
    }

    // Floor :
    push_object(new_Plane( 
        Plane(vec3(0,0,0),vec3(2,0,0),vec3(0,0,2)),
        COLOR_GROUND,0.0
    ));

}

////////////////////////////////////////////////////////////////////////
//                                                                    //
//                           Functions.glsl                           //
//                                                                    //
//                 Defining some useful math functions                //
//                                                                    //
////////////////////////////////////////////////////////////////////////

// 3d Matrix of rotation
mat3 rot(float t12, float t13, float t23){
    mat3 r1 = mat3(
        cos(t12),sin(t12),0,
        -sin(t12),cos(t12), 0.,
        0., 0., 1.
    );
    
    mat3 r2 = mat3(
        cos(t13), 0., sin(t13),
        0., 1., 0.,
        -sin(t13), 0., cos(t13)
    );
    
    mat3 r3 = mat3(
        1., 0., 0.,
        0., cos(t23), sin(t23),
        0., -sin(t23), cos(t23)
    );
    
    
    return r1 * r2 * r3;
}

// Norm squared
float norm2(vec3 v){
    return dot(v,v);
}

float norm(vec3 v){
    return sqrt(dot(v,v));
}

vec3 normalize(vec3 v){
    return v / norm(v);
}

float dist(vec3 A, vec3 B){
    return sqrt(dot(A-B,A-B));
}

// Distance squared
float dist2(vec3 A, vec3 B){
    return dot(A-B,A-B);
}

// Projection of the point A on the line L
vec3 projection(Line L, vec3 A){
    return L.origin + dot(A-L.origin,L.v) * L.v / norm2(L.v) ;
}

// Smooth minimum of a and b
float smooth_min(float a, float b){
    // return a < b ? a : b ;
    float k = 0.02 ;
    float h = a-b;
    return 0.5*( (a+b) - sqrt(h*h+k) );
}


// A pseudo-random number between 0 and 1
float random(float x){
    return fract(sin(x)*10000.0) ;
}

// A pseudo-random vec3
vec3 random3(vec3 v){
    return 2*vec3(
        //random(x*34.4938+836.9372),
        //random(x*02.3847+972.3085),
        //random(x*82.2984+184.3234)
        random(dot(v,vec3(0.1,294.284,-5.4824))+19.9382),
        random(dot(v,vec3(192.2482,-28.249,94.24))-21.2082),
        random(dot(v,vec3(-2.192,13.282,20.8294))+294.2972)
    ) - vec3(1) ;
}



// Get the decomposition of v in the basis (u1,u2)
vec3 locals_cord(vec3 v, vec3 u1, vec3 u2){
    // return cords of v in local base (u1,u2)
    float a = dot(u2,u2) * dot(v,u1) - dot(u1,u2) * dot(v,u2);
    float b = dot(u1,u1) * dot(v,u2) - dot(u1,u2) * dot(v,u1);
    float k = dot(u1,u1) * dot(u2,u2) - dot(u1,u2) * dot(u2,u1);
    return vec3(a,b,k) ;
}

// Return the spherical coordinate for a normed vector
vec2 spherical_cord(vec3 u){
    vec3 v = normalize(u) ;
    float angle = atan(v.z,v.x) / TWO_PI ;
    if (angle < 0) angle += 1 ;

    return vec2(angle, (1-v.y)/2 ) ;
    // angles in the range [-1,1]
}

////////////////////////////////////////////////////////////////////////
//                                                                    //
//                            Optical.glsl                            //
//                                                                    //
//                      Implementing Snell's laws                     //
//                                                                    //
////////////////////////////////////////////////////////////////////////

// Return the reflexion of v with respect to the normal n
vec3 reflexion(vec3 v, vec3 n){
    return v - 2.0 * dot(v,n) * n / dot(n,n) ;
}
// 
vec3 refraction(vec3 v_, vec3 n, float r){
    ///////////////////
    //  r = n1 / n2  //
    ///////////////////
    vec3 v = normalize(v_);

    // float n1 = r ;
    // float n2 = 1 ;

    float a = dot(n,n) ;
    float b = r * dot(n,v) ;
    float c = r * r * dot(v,v) - 1 ;

    if (b * b - a * c < 0) return vec3(0.);
    float d = sqrt(b * b - a * c);
    if ( dot(v,n) < 0) d*= -1 ;
    float k = (- b + d ) / a ;

    return k * n + r * v ;

    // float lambda = 1 - n1*n1 / (n2*n2) * dot(v,n) * dot(v,n) ;
    // return -r*n + sqrt(lambda) * v ;
}

vec2 taux_ref_old(float i1, float i2,float r){
    // https://claude-gimenes.fr/physique/propagation-des-ondes-electro-magnetiques/-ii-reflexion-et-refraction-en-milieux-isotropes-1
    // https://fr.wikipedia.org/wiki/Coefficient_de_Fresnel
    // http://www.joelsornette.fr/Archives/exotypes/exotype26.pdf

    // float r = (i1 - i2) / (i1 + i2) ;
    // float t = 2 * i1 / (i1 + i2) ;

    // float r = tan(i2-i1) / tan(i1+i2) ;
    // float t = 2 * sin(i2) * cos(i1) / (sin(i1+i2) * cos(i2-i1)) ;

    // float r = (cos(i1) - r*cos(i2)) / (cos(i1) + r*cos(i2));
    // float t = 2 * cos(i1) / (cos(i1) + r*cos(i2));

    // if (r < 1 && sqrt(1-r*r) >= 1) {return vec2(1,0);}

    // float n1 = 1 ;
    // float n2 = r * n1 ;

    // float R_ = (n2 * cos(i1) - n1 * cos(i2)) / (n2 * cos(i1) + n1 * cos(i2)) ;
    // float R = R_ * R_ ;
    // float T_1 = 4 * n1 * n2 * cos(i1) * cos(i2);
    // float T_2 = n2 * cos(i1) + n1 * cos(i2) ;

    // float T = T_1 / (T_2 * T_2);

    // return vec2(R,T);
    return vec2(0,0);

    // R(x) = (n2 * cos(x) - n1 * cos(I(x))) ** 2 / (n2 * cos(x) + n1 * cos(I(x))) ** 2
    // T(x) = 4 * n1 * n2 * cos(x) * cos(I(x)) / (n2 * cos(x) + n1 * cos(I(x))) ** 2 
}

vec2 taux_ref(vec3 u, vec3 v_t, vec3 n, float r){
    // http://www.joelsornette.fr/Archives/exotypes/exotype26.pdf
    // u Light vector, v_t transmitted vector, n normal vector, r quotient of refraction

    // if (r < 1 && sqrt(1-r*r) >= 1){
    //     return vec2(0,1);
    // }

    float n1 = 1 ;
    float n2 = r * n1 ;

    float i1 = dot(u,n) ;
    float i2 = dot(v_t,-n) ;

    float R_1 = n2 * i1 - n1 * i2 ;
    float R_2 = n2 * i1 + n1 * i2 ;
    float R_3 = R_1 / R_2 ;
    float R = R_3 * R_3 ;

    if (R_2 == 0.0) return vec2(1,0);

    float T_1 = 4 * n1 * n2 * i1 * i2;
    float T_2 = n2 * i1 + n1 * i2 ;
    float T = T_1 / T_2 ;

    return vec2(T,R);
}

////////////////////////////////////////////////////////////////////////
//                                                                    //
//                              SDF.glsl                              //
//                                                                    //
//                   Defining SDF for all primitives                  //
//                                                                    //
////////////////////////////////////////////////////////////////////////

float MandelbulbSDF(vec3 pos, float pow_, const int max_itr) {
    // https://editor.p5js.org/Taxen99/sketches/47CDg5-nV

    vec3 zeta = pos; 
    float dr = 1.0; // magic variable
    float r = 0.0; // the radius
    float theta ;
    float phi ;

    for (int n = 0; n > -1; n++) {
        if (n > max_itr) break;
        // vec3 v = zeta ;
        r = norm(zeta) ;
        theta = atan(sqrt(dot(zeta.xy,zeta.xy)), zeta.z);
        phi = atan(zeta.y, zeta.x);

        if (r > 2.0) break;

        dr = pow(r, pow_ - 1.0) * pow_ * dr + 1.0; // magic formula

        vec3 powered ;
        powered.x = pow(r, pow_) * sin(theta * pow_) * cos(phi * pow_);
        powered.y = pow(r, pow_) * sin(theta * pow_) * sin(phi * pow_);
        powered.z = pow(r, pow_) * cos(theta * pow_);
        // raise everything to the power of pow_

        zeta = powered + pos;
    }

    return 0.5 * log(r) * r / dr; // more magic to compute distance
}


float SDF(vec3 P, Sphere S){
    return dist(P, S.center) - S.radius ;
}

float SDF(vec3 P, Line L){
    vec3 P2 = projection(L,P);
    float a = dist(P2,L.origin) ;
    if (a < 0) return dist(L.origin,P) ;
    if (a > norm(L.v)) return dist(L.origin+L.v,P) ;
    return dist(P,P2) ;
}

float SDF(vec3 P, Cylinder C){
    vec3 P2 = projection(Line(C.origin,C.v), P) ;

    float h = dist(P2, C.origin) ;
    float r = dist(P2, P) ;
    
    float h2 = clamp(h,0,norm(C.v)) ;
    float r2 = C.radius ; // ??????

    vec3 P3 = C.origin + h2 * normalize(C.v) + r2 * normalize(P-P2) ;
    return dist(P,P3) ; 
}

float SDF(vec3 P, Plane T, int type){
    vec3 pos = locals_cord(P-T.origin,T.u1,T.u2) ;
    float a = pos.x / pos.z ;
    float b = pos.y / pos.z ;

    float a2, b2 ;

    if (type == TYPE_PLANE){
        a2 = clamp(a,-100,100) ;
        b2 = clamp(b,-100,100) ;
    } else if (type == TYPE_RECTANGLE){
        a2 = clamp(a,0,1) ;
        b2 = clamp(b,0,1) ;
    } else if (type == TYPE_CIRCLE){
        float r = a*a+b*b;
        if (r > 1){
            a2 = a / sqrt(r) ;
            b2 = b / sqrt(r) ;
        } else {
            a2 = a ;
            b2 = b ;
        }
    } else if (type == TYPE_TRIANGLE){
        float d = a+b ;
        if (d>1){
            a2 = a / d;
            b2 = b / d;
        } else {
            a2 = a ;
            b2 = b ;
        }
    } else {
        return -1. ;
    }

    vec3 P2 = T.origin + a2 * T.u1 + b2 * T.u2 ;
    return dist(P,P2) ;
}

float SDF(Line L, Object O){
    vec3 P = L.origin ;
    if (O.type == TYPE_ERROR) return -1.;
    if (O.type == TYPE_SPHERE) return SDF(P, O.sphere);
    if (O.type == TYPE_CYLINDER) return SDF(P, O.cylinder);
    if (O.type == TYPE_MANDELBULB) return MandelbulbSDF(P-O.sphere.center, 5.0,30);
    return SDF(P, O.plane, O.type);
}

float SDF(Line L, Group G){
    float min_dist = -1. ;
    for (int j = 0 ; j < G.nb_objects ; j++){
        float dist_j = SDF(L, G.objects[j]) ;
        if ( dist_j < 0 ) continue ;
        if (min_dist == -1.0){
            min_dist = dist_j ;
        } else {
            if (G.type == STRUCT_UNION) min_dist = min(min_dist, dist_j) ;
            if (G.type == STRUCT_SMOOTH_U) min_dist = smooth_min(min_dist, dist_j) ;
            if (G.type == STRUCT_INTERSECTION) min_dist = max(min_dist, dist_j) ;
        }
    }
    return min_dist ;
}

Object SDF(Line Ray){
    float min_dist = -1. ;
    Object best_obj ;
    

    for (int j = 0 ; j < nb_object ; j++){
        float dist_j = SDF(Ray, objects[j]) ;
        if ( dist_j < 0 ) continue ;
        if (min_dist == -1.){
            min_dist = dist_j ;
            best_obj = objects[j] ;
        } else {
            best_obj = (dist_j < min_dist) ? objects[j] : best_obj ;
            min_dist = smooth_min(min_dist, dist_j) ;
        }
    }
    return min_dist == -1 ? null_Object() : best_obj ;
}

////////////////////////////////////////////////////////////////////////
//                                                                    //
//                         Intersections.glsl                         //
//                                                                    //
//                   Calculating the intersections                    //
//                 between a light ray and an object                  //
//                                                                    //
////////////////////////////////////////////////////////////////////////


Line get_intersection(Line L, Sphere S){
    // Return the intersection and the normal of the sphere
    vec3 Lv = normalize(L.v) ;
    float d = dot(S.center - L.origin,Lv) ;
    if (d <= 0.0) return Line(vec3(0.),vec3(0.));
    vec3 H = L.origin + d * Lv;
    vec3 u = H - S.center;
    
    float l = S.radius * S.radius - dot(u,u);

    // H isnt in the sphere
    if (l < 0.) return Line(H,vec3(0.));
    
    // H is in the sphere
    // We calculate the intersection based on H
    
    vec3 I = H ;
    if (d - sqrt(l) < 0) I += sqrt(l) * Lv;
    else I -= sqrt(l) * Lv;
    // if (dot(I-L.origin, L.v) < 0) I = H - sqrt(l) * Lv;

    vec3 n = I - S.center;

    return Line(I, n) ;
}

Line get_intersection(Line L, Cylinder C){

    vec3 v1 = cross(L.v,C.v);
    vec3 v2 = cross(L.origin - C.origin, C.v);

    
    float a = norm2(v1) ;
    float b = dot( v1, v2 );
    float c = norm2( v2 ) - C.radius * C.radius * norm2(C.v) ;

    float d = b*b - a*c ;
    float lambda ;

    if (d < 0.) return Line(vec3(0.),vec3(0.));
    else if (d == 0.0) lambda = - b / a ;
    else {
        float lambda1 = (- b + sqrt(d)) / a ;
        float lambda2 = (- b - sqrt(d)) / a ;

        if ( lambda1 < 0 && lambda2 < 0) return Line(vec3(0.),vec3(0.));
        else if (lambda1 < 0) lambda = lambda2;
        else if (lambda2 < 0) lambda = lambda1;
        else if (lambda1 < lambda2) lambda = lambda1;
        else lambda = lambda2;
    }
    
    vec3 I = L.origin + lambda * L.v ;
    vec3 N = I - C.origin - dot(I - C.origin, C.v) * C.v / norm2(C.v) ;

    if (
        dot(I-C.origin, C.v) < 0.0 
        || dot(I-C.origin, C.v) > norm2(C.v)
        || dot(I-L.origin, L.v) < 0.0
    ) return Line(vec3(0.),vec3(0.));

    return Line(I,N);
}

Line get_intersection(Line L, Plane T, int type){

    vec3 n = cross(T.u1,T.u2);
    
    // Line is perpendicular to n --> no intersection
    if (dot(n,L.v) == 0.0) return Line(vec3(0.),vec3(0.));

    // The line intersect the plane of the plane
    float lambda = dot(n, T.origin - L.origin) / dot(n,L.v);
    if (lambda <= 0.0) return Line(vec3(0.),vec3(0.)); // The intersection is behind us
    
    // H is the intersection
    vec3 H = L.origin + lambda * L.v ;
    vec3 v = H - T.origin ;

    vec3 normal = dot(L.v,n) > 0.0 ? -n : n ;
    normal = normalize(normal) ;

    if (type == TYPE_PLANE) return Line(H, normal);

    // a and b are coordinate relative to the local base of the plane (T.u1, T.u2)
    vec3 coordinate = locals_cord(v,T.u1,T.u2);
    if (coordinate.z == 0.0) return Line(vec3(0.),vec3(0.));;
    float a = coordinate.x / coordinate.z;
    float b = coordinate.y / coordinate.z; 
    

    // The intersection is in the rectangle
    if (type == TYPE_RECTANGLE){
        if (0. <= a && a <= 1. && 0. <= b && b <= 1.){
            return Line(H, normal);
        } else {
            return Line(vec3(0.),vec3(0.));
        }
    }

    // The intersection is in the circle
    if (type == TYPE_CIRCLE){
        if ( a*a + b*b < 1){
            return Line(H, normal);
        } else {
            return Line(vec3(0.),vec3(0.));
        }
    }



    // The intersection is in the triangle
    if (0. <= a+b && a+b <= 1. && 0. <= a && 0. <= b){
        return Line(H, normal);
    } else { return Line(vec3(0.),vec3(0.));}

}

Line get_intersection(Line L, Object O){
    if (O.type == TYPE_ERROR) return Line(vec3(0.),vec3(0.));
    if (O.type == TYPE_SPHERE) return get_intersection(L, O.sphere);
    if (O.type == TYPE_CYLINDER) return get_intersection(L, O.cylinder);
    if (O.type == TYPE_MANDELBULB) return get_intersection(L, Sphere(O.sphere.center,.3));
    return get_intersection(L, O.plane, O.type);
}

Object get_intersection(Line Ray){
    float min_len = -1. ;
    Object best_obj ;

    for (int j = 0 ; j < nb_object ; j++){
        Line Intersection = get_intersection(Ray,objects[j]);
        float dist = dot(Intersection.origin - Ray.origin,Intersection.origin - Ray.origin);
        if (
            Intersection.v != vec3(0.) && 
            abs(dist) > 0.001  &&
            ( min_len == -1. || dist < min_len )
        ){
            min_len = dist ;
            best_obj = objects[j] ;
        }
    }
    return min_len == -1 ? null_Object() : best_obj ;
}


////////////////////////////////////////////////////////////////////////
//                                                                    //
//                             Color.glsl                             //
//                                                                    //
//                  Calculating the color of a point                  //
//                                                                    //
////////////////////////////////////////////////////////////////////////


vec3 img_color(sampler2D img, Object Obj, Line Normal){

    vec2 pos ;
    if (Obj.type == TYPE_SPHERE){
        pos = spherical_cord(Normal.origin - Obj.sphere.center);
    } else if (Obj.type == TYPE_CYLINDER){
        vec3 v = Normal.origin - Obj.cylinder.origin ;
        vec2 pos = vec2(0,0);
        pos.y = .5 ;
        // pos.y = dot(v,Obj.cylinder.v) / sqrt(norm2(Obj.cylinder.v)) ;
        // v -= Obj.cylinder.radius * pos.y * Obj.cylinder.v ;
        pos.x = atan(v.z,v.x) / TWO_PI ;
    } else {
        vec3 v = Normal.origin - Obj.plane.origin ;
        vec3 v2 = locals_cord(v,Obj.plane.u1,Obj.plane.u2);
        pos = vec2(v2.x / v2.z, 1-v2.y / v2.z) ;
    }

    // return pos ;
    return texture2D(img,pos).rgb ;
}

vec3 lattice_color(vec2 pos, vec3 c1, vec3 c2) {
    bool cond = (mod(pos.x,2) > 1) == (mod(pos.y,2) > 2);
    return cond ? c1 : c2 ;
}

vec4 calc_color(Object Obj, Line Normal){

    vec3 col = vec3(0);
    vec3 obj_col ;
    float n_col = 0.0;

    if        ( Obj.color == COLOR_HOKUSAI){
        obj_col = img_color(hokusai, Obj, Normal);
        // obj_col = texture2D(hokusai,img_color(Obj, Normal)).rgb ;

    } else if ( Obj.color == COLOR_GALAXY){
        obj_col = img_color(galaxy, Obj, Normal);
    } else if ( Obj.color == COLOR_EARTH){
        obj_col = img_color(earth, Obj, Normal);
    } else if ( Obj.color == COLOR_SUN){
        obj_col = img_color(sun, Obj, Normal);

    } else if ( Obj.color == COLOR_MOON){
        obj_col = img_color(moon, Obj, Normal);
    } else if ( Obj.color == COLOR_SAGITTARIUS){
        obj_col = img_color(sagittarius_A, Obj, Normal);

    } else if ( Obj.color == COLOR_GROUND) {
        vec2 pos = fract(Normal.origin.xz / 10) ;
        if (draw_type == 2 && false) obj_col = vec3(1,0,0) ;
        else obj_col = texture2D(ground,pos).rgb ;
    } else if ( Obj.color == COLOR_LATTICE_1) {
        obj_col = lattice_color(Normal.origin.xz,vec3(0.25),vec3(0.75)) ;
    } else if ( Obj.color == COLOR_LATTICE_2) {
        obj_col = lattice_color(Normal.origin.xz,vec3(0),vec3(1)) ;
    } else if ( Obj.color == COLOR_LATTICE_3) {
        obj_col = lattice_color(Normal.origin.xz,vec3(1,0,0),vec3(1, 1, 0)) ;
    } else if ( Obj.color == COLOR_LATTICE_4) {
        obj_col = lattice_color(Normal.origin.xz,vec3(0,0,1),vec3(0,1,0)) ;
    } else if ( Obj.color == COLOR_LATTICE_5) {
        obj_col = lattice_color(Normal.origin.xz,vec3(0,1,1),vec3(0, 1, 0.5)) ;

    } else { // base color
        obj_col = Obj.color ;
        // obj_col = normalize(Normal.v) ;
    }


    if (Obj.is_light){
        return vec4(obj_col,-1.0);
    } else if (light_type != 0 && Obj.mirror != 1.0) { // Add light
        for (int i=0 ; i < nb_object ; i++) {

            Object Light = objects[i];
            if ( ! Light.is_light ){ continue; }

            vec3 center ; 
            if (Light.type == TYPE_SPHERE) center = Light.sphere.center ; 
            else if (Light.type == TYPE_CYLINDER) center = Light.cylinder.origin ;
            else center = Light.plane.origin ;

            vec3 light_dir = center - Normal.origin ;
            Line ray = Line(center,-light_dir);
            Object inter = get_intersection(ray);

            if (
                inter.type == Obj.type &&
                inter.color == Obj.color &&
                inter.mirror == Obj.mirror &&
                inter.is_light == Obj.is_light &&
                inter.sphere == Obj.sphere &&
                inter.plane == Obj.plane
            ){
                // angle : 1 if normal    0 if colinear
                float angle = dot(normalize(Normal.v),normalize(light_dir));
                angle = clamp(angle,0.,1.);
                float dist = sqrt(norm2(light_dir));
                // float dist = norm2(light_dir);
                float coef ;
                coef = angle * 5 / (2.0 + dist) ;
                
                if (angle != 0){
                    vec3 Light_color = vec3(1);
                    // vec3 Light_color = Light.color;
                    if (light_type == 2) col += coef * Light_color * obj_col;
                    else col += coef * Light_color;
                    n_col += coef ;
                }
            } else {
                col += vec3(0);
                n_col += 1 ;
            }

        }
    }

    if (light_type != 2) { // Ambient light
        // col += ( 1 - Obj.mirror ) * obj_col ;
        // n_col += ( 1 - Obj.mirror ) ;
        col += obj_col ;
        n_col += 1 ;
    }
    
    return vec4(col,n_col);
}

////////////////////////////////////////////////////////////////////////
//                                                                    //
//                            Rendering.glsl                          //
//                                                                    //
//                          Render the scene                          //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#define NB_ITER 6
#define NB_ITER_3 15

vec3 draw(Line Ray){
    vec4 col = vec4(0.);
    float col_weight = 1;

    for (int i=0 ; i<NB_ITER ; i++){

        Object best_obj = get_intersection(Ray);

        if (best_obj.type == TYPE_ERROR){
            if (light_type == 2) col += col_weight * vec4(0,0,0,1) ;
            else col += col_weight * vec4(0.0, 0.0, 0.3,1) ;
            break ;
        }

        Line N = get_intersection(Ray,best_obj);

        vec4 res = calc_color(best_obj, N );
        
        if ( res.w == -1.0) { // The object is light
            col += col_weight * vec4(res.xyz,1);
            break ;
        }

        col += col_weight * res.w * vec4(res.xyz,1) ;
        col_weight *= best_obj.mirror ;

        if (col_weight == 0) break ;
        
        if (best_obj.refrac != 1.0 ){
            bool is_reverse ;
            if (best_obj.type == 0) is_reverse = dot(N.origin - best_obj.sphere.center,Ray.v) > 0.0 ;
            else if (best_obj.type == 1) is_reverse = dot(N.v,Ray.v) > 0.0 ;
            else is_reverse = dot(cross(best_obj.plane.u1,best_obj.plane.u2),Ray.v) > 0 ;
            
            float r = is_reverse ? 1/best_obj.refrac : best_obj.refrac ;

            Ray.v = refraction(Ray.v,N.v,r) ;
            if (Ray.v == vec3(0.0)) break ;
        } else {
            break ;
            Ray.v = reflexion(Ray.v,N.v);
        }
        Ray.origin = N.origin ;

        // iter loop for calc reflexion or refraction
    }

    return col.rgb / col.w ;
}

vec3 draw2(Line Ray){

    Line Rays[NB_ITER] ;
    vec3 cols[NB_ITER] ;
    float coeffs[NB_ITER] ;
    bool is_init[NB_ITER] ;

    // Init Array
    for (int i = 0; i < NB_ITER; i++){
        Rays[i] = Line(vec3(0),vec3(0)) ;
        cols[i] = vec3(0) ;
        coeffs[i] = 1 ;
        is_init[i] = false ;
    }

    // Init first Ray
    int n = 0 ;
    Rays[n] = Ray ;
    is_init[n] = true ;
    n++ ;


    for (int i = 0; i < NB_ITER; i++){

        // if (i>light_height) break ;
        if (! is_init[i]) break ;
        Object best_obj = get_intersection(Rays[i]);

        if (best_obj.type == TYPE_ERROR){ // show sky
            if (light_type == 2) cols[i] =  vec3(0.0) ;
            else {
                vec2 v = spherical_cord(Rays[i].v + vec3(0,.1,0)) ;
                if (v.y > .5) cols[i] = vec3(0.0, 0.8, 1.0) ;
                else cols[i] = texture2D(sky,vec2(v.x,2*v.y)).rgb ;
            }
            continue ;
        }

        Line N = get_intersection(Rays[i],best_obj);

        vec4 res = calc_color(best_obj, N );
        
        if ( res.w == -1.0) { // obj is light
            cols[i] = res.rgb;
            continue ;
        }
        cols[i] = res.rgb ;
        // *0 + vec3(0.0, 0.5176, 1.0)
        if ( n-2< NB_ITER && best_obj.refrac != 1){

            bool is_reverse ;
            if (best_obj.type == TYPE_SPHERE)
                is_reverse = dot(N.origin - best_obj.sphere.center,Ray.v) > 0.0 ;
            else if (best_obj.type == TYPE_CYLINDER)
                is_reverse = dot(N.v,Ray.v) > 0.0 ;
            // else if (best_obj.type == 1) is_reverse = dot(N.origin - best_obj.cylinder.origin,Ray.v) > 0.0 ;
            else is_reverse = dot(cross(best_obj.plane.u1,best_obj.plane.u2),Ray.v) > 0 ;
            float r = is_reverse ? 1/best_obj.refrac : best_obj.refrac ;
            r = abs(r) ;

            Rays[n].v = refraction(Rays[i].v,N.v,r);        // Refraction
            Rays[n].origin = N.origin ;
            vec2 C = taux_ref(Rays[i].v,Rays[n].v,N.v,r) ;
            coeffs[n] = C.x * coeffs[i]/res.w;
            is_init[n] = true ;
            n++ ;

            Rays[n].v = reflexion(Rays[i].v,N.v);          // Reflexion
            Rays[n].origin = N.origin ;
            coeffs[n] = C.y * coeffs[i]/res.w;
            is_init[n] = true ;
            n++ ;

        } else if (best_obj.mirror != 0.0 && n<NB_ITER) {
            Rays[n].v = reflexion(Rays[i].v,N.v);
            Rays[n].v += random3(Rays[i].v) * 0.0;
            Rays[n].origin = N.origin ;
            coeffs[n] = coeffs[i] * best_obj.mirror;
            is_init[n] = true ;
            n++ ;
        }
        coeffs[i] *= (1-best_obj.mirror) * res.w ;
    }

    vec3 col = vec3(0) ;
    float sum_coeffs = 0 ;
    int nb_init = 0 ;

    for (int i = 0; i < NB_ITER; i++){
        if (is_init[i]){
            col += cols[i] * coeffs[i] ; 
            sum_coeffs += coeffs[i] ;
            // sum_coeffs += 1 ;
            nb_init += 1 ;
        }
    }
    
    if (sum_coeffs == 0) return vec3(0.8039, 0.149, 0.7686) ;
    return col / sum_coeffs ;
}

vec3 draw3(Line Ray){
    vec4 col = vec4(0);
    
    float dist = 0.0 ;
    
    Object best_obj;
    vec3 obj_col;

    for (int i = 0 ; i < NB_ITER_3 ; i++){
        best_obj = SDF(Ray);

        float D = (best_obj.type == TYPE_ERROR) ? -1 : SDF(Ray,best_obj); 

        obj_col = calc_color(best_obj,Line(Ray.origin,vec3(0))).rgb ;
        // obj_col = best_obj.color ;

        if (abs(D) < .1){
            // float dist_max = 30 ;
            // float alpha_max = 0.5 ;
            // float alpha = dist > dist_max ? 1 : (1 - pow(1-dist/dist_max,2)) ;
            // col += vec4(0.0, 0.0, 0.0, alpha_max*alpha) ;
            col += vec4(obj_col,1) ;
            break ;
            // return col.rgb / col.w ;
        }
        // col += vec4(obj_col,1) * exp(-D/10) ;
        // col += vec4(obj_col,1) * exp(-D*D/25) ;
        // col += vec4(1.0, 1.0, 1.0, 1.0) * exp(-D*D/40) ;

        if (D < 0) return vec3(.5,.5,.5) ; // Inside something
        if (D > 20){ // show sky
            col += vec4(0.0, 0.0, 0.3,10);
            return vec3(0.0, 0.0, 0.3);
            return col.rgb / col.w ;
            // if (light_type == 2) cols[i] =  vec3(0.0) ;
            // else {
            //     vec2 v = spherical_cord(Rays[i].v + vec3(0,.1,0)) ;
            //     if (v.y > .5) cols[i] = vec3(0.0, 0.8, 1.0) ;
            //     else cols[i] = texture2D(sky,vec2(v.x,2*v.y)).rgb ;
            // }
        }

        Ray.origin += D * Ray.v ;
        dist += D ;
    }
    
    // Do I intersect something ?
    vec3 base_color = vec3(1.0, 0.0, 0.9) ;
    //vec3 base_color = obj_col ;
    return (col.w==0) ? base_color : col.rgb / col.w;
}

////////////////////////////////////////////////////////////////////////
//                                                                    //
//                              main.glsl                             //
//                                                                    //
//                     Caculate the position and                      //
//                  orientation to render the scene                   //
//                                                                    //
////////////////////////////////////////////////////////////////////////

void main() {
    float size = max(u_resolution.x,u_resolution.y);
    vec2 st =  1.0 * (gl_FragCoord.xy / vec2(size,size) - .5);

    init_objects();

    vec2 mouse = 3 * (u_mouse / u_resolution - .5) ;
    mat3 M = rot(0,-2*mouse.x,mouse.y);
    vec3 B = vec3(st.x,st.y,0.66);
    vec3 A = vec3(position.x,2,position.y);

    Line Ray = Line(
        A,
        normalize(M * B)
    );
    vec3 col ;
    if      (draw_type == 0) col = draw(Ray) ;
    else if (draw_type == 1) col = draw2(Ray) ;
    else if (draw_type == 2) col = draw3(Ray);
    else col = vec3(0.8039, 0.149, 0.7686) ;
    
    gl_FragColor = vec4(col, 1.0);
   
}


