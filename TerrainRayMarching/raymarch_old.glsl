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

// Terrain or Ray Marching :
// https://blog.maximeheckel.com/posts/painting-with-math-a-gentle-study-of-raymarching/
// https://www.shadertoy.com/view/XlyXRh Simple Ray Marching
// https://www.shadertoy.com/view/4ttSWf Rainforest
// https://www.shadertoy.com/view/3dGSWR FBM Terrain

// https://www.youtube.com/watch?v=BFld4EBO2RE Video Painting a Landscape with Maths (IQ)


// #ifdef GL_ES
// precision mediump float;
precision highp float;
// #endif

// uniform sampler2D sagittarius_A;

uniform bool bases_vectors ;
uniform int scene ;

uniform vec2 u_resolution;
uniform vec2 u_mouse;
uniform float u_time;
uniform vec2 position;
uniform int variable;


#define PI 3.1415926535897932384626433832795
#define TWO_PI 6.283185307179586476925286766559

#define COLOR_LATTICE_1 vec3(0,0.01,0)
#define COLOR_LATTICE_2 vec3(0,0.02,0)
#define COLOR_LATTICE_3 vec3(0,0.03,0)
#define COLOR_LATTICE_4 vec3(0,0.04,0)
#define COLOR_LATTICE_5 vec3(0,0.05,0)


#define TYPE_ERROR -1
#define TYPE_SPHERE 1
#define TYPE_CYLINDER 2
#define TYPE_PLANE 3
#define TYPE_RECTANGLE 4
#define TYPE_CIRCLE 5
#define TYPE_TRIANGLE 6
#define TYPE_MANDELBULB 7
#define TYPE_TERRAIN 8


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


////////////////////////////////////////////////////////////////////////

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

Object new_Object(Sphere S, vec3 color, float mirror){
    Object O = null_Object();
    O.sphere = S;
    O.type = TYPE_SPHERE;
    O.color = color;
    O.mirror = mirror;
    return O;
}

Object new_Object(Cylinder C, vec3 color, float mirror){
    Object O = null_Object();
    O.cylinder = C;
    O.type = TYPE_CYLINDER;
    O.color = color;
    O.mirror = mirror;
    return O;
}

Object new_Object(Plane P, vec3 color, float mirror){
    Object O = null_Object();
    O.plane = P;
    O.type = TYPE_PLANE;
    O.color = color;
    O.mirror = mirror;
    return O;
}
Object new_Rectangle(Plane P, vec3 color, float mirror){
    Object O = new_Object(P, color, mirror);
    O.type = TYPE_RECTANGLE;
    return O;
}
Object new_Plane(Plane P, vec3 color, float mirror){
    Object O = new_Object(P, color, mirror);
    O.type = TYPE_PLANE;
    return O;
}
Object new_Circle(Plane P, vec3 color, float mirror){
    Object O = new_Object(P, color, mirror);
    O.type = TYPE_CIRCLE;
    return O;
}
Object new_typed_Plan(Plane P, vec3 color, float mirror, int type){
    Object O = new_Object(P, color, mirror);
    O.type = type;
    return O;
}

Object new_Bulb(vec3 pos, vec3 color, float mirror){
    Object O = null_Object();
    O.sphere = Sphere(pos,0);
    O.type = TYPE_MANDELBULB;
    O.color = color;
    O.mirror = mirror;
    return O;
}

Object new_terrain(){
    Object O = null_Object();
    O.type = TYPE_TERRAIN;
    O.color = vec3(0.0, 0.5, 0.0);
    return O;
}

Object Refractor(Object O, float refrac){
    O.refrac = refrac;
    return O;
}

Object Light(Object O){
    O.is_light = true;
    return O ;
}



////////////////////////////////////////////////////////////////////////

#define nb_object 15

// Idea :
// List of objects
// List of links
// obj : [0,1,2,3,4,5]
// links : [union 0 1, inter 2 3 , ...]
// Calc : [o0, o1, o2, o3, o4, o5, union 0 1, inter 7 2, ...]
// Calc : Array[3*n]

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

void push_object(Object obj){
    // inout int index_obj, 
    objects[index_obj] = obj ;
    index_obj++ ;
}

void push_cube(vec3 pos,vec3 u1,vec3 u2,vec3 u3, vec3 col, float mirror){
    push_object(new_Rectangle(Plane(pos,u1,u2),col,mirror));
    push_object(new_Rectangle(Plane(pos,u2,u3),col,mirror));
    push_object(new_Rectangle(Plane(pos,u3,u1),col,mirror));
    push_object(new_Rectangle(Plane(pos+u1+u2+u3,-u1,-u2),col,mirror));
    push_object(new_Rectangle(Plane(pos+u1+u2+u3,-u2,-u3),col,mirror));
    push_object(new_Rectangle(Plane(pos+u1+u2+u3,-u3,-u1),col,mirror));
}

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
            Sphere(vec3(0,variable/2.5,0),1.2),
            vec3(1.0, 0.8157, 0.0),
            0.0
        )));
        push_object(new_Rectangle(
            Plane(vec3(-5,1,5),vec3(10,0,0),vec3(0,6,0)),
            vec3(0.0, 0.0, 1.0), 0.0
        ));
        // push_object(new_Object(Sphere(vec3(-2,1.5,0),1),vec3(0.0, 1.0, 0.702),0.0));
        // push_object(new_Object(Sphere(vec3(2.5,3.5,-1),.4),vec3(0.6588, 0.0706, 0.6392),0.0));

        push_object(new_Object(Sphere(vec3(-2,1.5,0),1),vec3(0.5608, 0.4471, 0.2706),0.0));
        push_object(new_Object(Sphere(vec3(2.5,3.5,-1),.4),vec3(.5,.8,.4),0.0));
        
        // push_object(new_Object(Cylinder(vec3(-4,1,3),vec3(2,3,-1),2),vec3(0.1, 0.1, 1.0),0));
        // push_object(new_Object(Cylinder(vec3(0,0,0),vec3(0,1.01,0),1),COLOR_SKY,0));
    }
    

    // Scene 2 : Mirrors
    if (scene == 2){
        push_object(Light(new_Object(Sphere(vec3(0,4,0),.2),vec3(1.0, 0.8157, 0.0),0.0)));
        push_object(new_Rectangle(
           Plane(vec3(-7.5,0,7),vec3(15,0,0),vec3(0,10,0)),
           vec3(0.0, 0.0, 0.0),.9
        ));
        push_object(new_Object(Sphere(vec3(-2,1.5,0),1),vec3(0.5608, 0.4471, 0.2706),0.0));
        push_object(new_Object(Sphere(vec3(2,2.2,1),1.8),vec3(1.0, 0.6, 0.0),.9));
        push_object(new_Object(Cylinder(vec3(-4,1,1),vec3(2,3,-1),.5),vec3(0.0, 0.7, 1.0),.9));
    }
    
    // Scene 3 : Refraction
    if (scene == 3){

        // float refrac = variable / 10.0 - 1. + 0.01 ;
        float refrac = exp(variable/5) + 0.001;

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

    // Scene 4 : A sphere
    if (scene == 4){
        push_object(new_Bulb(vec3(-1,1,0),vec3(1),0)) ;
        vec3 col = vec3(0.3333, 0.7176, 0.5294);
        push_object(new_Object(Sphere(vec3(2,variable,0),1.5),col,0.0));
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

    // Scene 6 : A room
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

    // Scene 0 : A terrain
    if (scene == 0){
        push_object(new_terrain());
        return ;
    }

    // Floor :
    push_object(new_Plane( 
        Plane(vec3(0,0,0),vec3(2,0,0),vec3(0,0,2)),
        vec3(0.4941, 0.3647, 0.1804),0.0
    ));

}

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

float dist2(vec3 A, vec3 B){
    return dot(A-B,A-B);
}

vec3 projection(Line L, vec3 A){
    return L.origin + dot(A-L.origin,L.v) * L.v / norm2(L.v) ;
}

float smooth_min(float a, float b){
    // return a < b ? a : b ;
    float k = 0.02 ;
    float h = a-b;
    return 0.5*( (a+b) - sqrt(h*h+k) );
}

float random(float x){
    return fract(sin(x)*10000.0) ;
}

vec3 random3(float x){
    return 2*vec3(
        random(x*34.4938+836.9372),
        random(x*02.3847+972.3085),
        random(x*82.2984+184.3234)
    ) - vec3(1) ;
}




vec3 locals_cord(vec3 v, vec3 u1, vec3 u2){
    // return cords of v in local base (u1,u2)
    float a = dot(u2,u2) * dot(v,u1) - dot(u1,u2) * dot(v,u2);
    float b = dot(u1,u1) * dot(v,u2) - dot(u1,u2) * dot(v,u1);
    float k = dot(u1,u1) * dot(u2,u2) - dot(u1,u2) * dot(u2,u1);
    return vec3(a,b,k) ;
}

vec2 spherical_cord(vec3 u){
    vec3 v = normalize(u) ;
    float angle = atan(v.z,v.x) / TWO_PI ;
    if (angle < 0) angle += 1 ;

    return vec2(angle, (1-v.y)/2 ) ;
}



// Return the reflexion of v with respect to the normal n
vec3 reflexion(vec3 v, vec3 n){
    return v - 2.0 * dot(v,n) * n / dot(n,n) ;
}

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

float SDFTerrain(vec3 P){
    vec2 pos = P.xz ;
    return P.y - 0.0 + sin(10*pos.x) * cos(5*pos.y) * 0.5 ;
}

float SDF(Line Ray, Object O){
    vec3 P = Ray.origin ;
    if (O.type == TYPE_ERROR) return -1.;
    if (O.type == TYPE_SPHERE) return SDF(P, O.sphere);
    if (O.type == TYPE_CYLINDER) return SDF(P, O.cylinder);
    if (O.type == TYPE_MANDELBULB) return MandelbulbSDF(P-O.sphere.center, 5.0,30);
    if (O.type == TYPE_TERRAIN) return SDFTerrain(Ray.origin);
    return SDF(P, O.plane, O.type);
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
            min_dist = (dist_j < min_dist) ? dist_j : min_dist ;
            // min_dist = smooth_min(min_dist, dist_j) ;
        }
    }
    return min_dist == -1 ? null_Object() : best_obj ;
}


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

    vec3 obj_col ;
    float n_col = 0.0;

    if ( Obj.color == COLOR_LATTICE_1) {
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

    return vec4(obj_col,n_col);
}

#define NB_ITER 100


vec3 draw(Line Ray){
    vec4 col = vec4(0);
    
    float dist = 0.0 ;

    for (int i = 0 ; i < NB_ITER ; i++){
        Object best_obj = SDF(Ray);

        float D = (best_obj.type == TYPE_ERROR) ? -1 : SDF(Ray,best_obj); 

        vec3 obj_col = calc_color(best_obj,Line(Ray.origin,vec3(0,2,0))).rgb ;
        // vec3 obj_col = best_obj.color ;

        if (abs(D) < .02){
            // float dist_max = 30 ;
            // float alpha_max = 0.5 ;
            // float alpha = dist > dist_max ? 1 : (1 - pow(1-dist/dist_max,2)) ;
            // col += vec4(0.0, 0.0, 0.0, alpha_max*alpha) ;
            col += vec4(obj_col,1) ;
            // break ;
            return obj_col ;
            // return col.rgb / col.w ;
        }
        // col += vec4(obj_col,1) * exp(-D/10) ;
        // col += vec4(obj_col,1) * exp(-D*D/25) ;
        // col += vec4(1.0, 1.0, 1.0, 1.0) * exp(-D*D/40) ;

        if (D < 0) return vec3(.5,.5,.5) ;
        if (D > 100){ // show sky
            col += vec4(0.0, 0.0, 0.3,10);
            return col.rgb ; // / col.w ;
            return vec3(0.0, 0.0, 0.3);
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
    return (col.w==0) ? vec3(1.0, 0.0, 0.9) : col.rgb ;
}

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
    vec3 col = draw(Ray) ;
    
    gl_FragColor = vec4(col, 1.0);
   
}


