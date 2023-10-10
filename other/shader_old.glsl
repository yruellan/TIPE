// Author:
// Title:

#ifdef GL_ES
// precision mediump float;
precision highp float;
#endif

uniform sampler2D u_tex0;
uniform vec2 u_tex0Resolution;

uniform vec2 u_resolution;
uniform vec2 u_mouse;
uniform float u_time;
uniform vec2 position;
uniform float light_height;
uniform int light_type; // 0: no light, 1: light+ambience, 2: only light

// Structure definitions
struct Line {
    vec3 point ;
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
    int type ;
    // -1 error, 0 sphere, 1 plane, 2 rect, 3+ triangle
    vec3 color ;
    float mirror ;
    float refrac ;
    bool is_light ;
};

////////////////////////////////////////////////////////////////////////

Object null_Object(){
    return Object(
        Sphere(vec3(0.),0.),
        Plane(vec3(0.),vec3(0.),vec3(0.)),
        -1,                     // type
        vec3(0.0),            // color
        0.0,                    // mirror
        1.0,                    // refrac
        false                   // is_light
    ) ;
}

Object new_Object(Sphere S, vec3 color, float mirror, bool is_light){
    return Object(S, Plane(vec3(0.),vec3(0.),vec3(0.)), 0, color, mirror, 1, is_light);
}
Object new_Object(Sphere S, vec3 color, float mirror){
    return Object(S, Plane(vec3(0.),vec3(0.),vec3(0.)), 0, color, mirror, 1,false);
}

Object new_Object(Plane T, vec3 color, float mirror, bool is_light){
    return Object(Sphere(vec3(0.),0.), T, 1, color, mirror, 1, is_light);
}
Object new_Object(Plane T, vec3 color, float mirror){
    return Object(Sphere(vec3(0.),0.), T, 1, color, mirror, 1, false);
}

Object new_Rectangle(Plane T, vec3 color, float r){
    return Object(Sphere(vec3(0.),0.), T, 2, color, 0.8, r, false);
}

////////////////////////////////////////////////////////////////////////

#define nb_object 7


Object objects[nb_object] = Object[nb_object](
    new_Object(Sphere(vec3(0,light_height,0),.5),vec3(0.9, 0.85, 1.0),0.0,true),
    new_Object(Sphere(vec3(1,0,-0.5),.5),vec3(0.0, 0.7, 0.0),0.0),
    new_Object(Sphere(vec3(1,1,0.5),.5),vec3(0.8824, 1.0, 0.0),0.0),
    new_Object(Sphere(vec3(-1,0,-0.5),.5),vec3(0.99, 0.0, 0.0),0.0),
    new_Object(
    Plane(vec3(0,-2,0),vec3(2,0,0),vec3(0,0,2)),
    vec3(0.5, 0.5, 0.5),0.0
    ),
    // new_Rectangle(Plane(vec3(-5,-1,0),vec3(0,0,2),vec3(0,2,0)),vec3(0,1,1),0.8),
    new_Rectangle(Plane(vec3(-6,-1,0),vec3(0,2,0),vec3(0,0,2)),vec3(1,0,1),0.8),
    // new_Object(
    // Plane(vec3(0,0,15),vec3(2,0,0),vec3(0,2,0)),
    // vec3(0.7, 0.0, 0.0),0.8
    // ),
    // new_Object(
    // Plane(vec3(15,0,0),vec3(0,0,2),vec3(0,2,0)),
    // vec3(0.0, 0.0, 0.7),0.8
    // ),
    Object(
        Sphere(vec3(0.),0.), Plane(vec3(5,0,0),vec3(-1,0,-3),vec3(0,2,0)),
        2,vec3(0,0,0.01), 0.0, 1.0, false
    )
);

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


float rand(vec2 co){
    return fract(sin(dot(co, vec2(12.9898, 78.233))) * 43758.5453);
}

vec3 random_dir(vec3 v){
    // float x = 2. * rand(v) - 1.;
    // float y = 2. * rand(vec3(v.y,v.z,v.x)) - 1.;
    // float z = 2. * rand(vec3(v.z,v.x,v.y)) - 1.;
    // return normalize(vec3(x,y,z));
    return vec3(0.);
}

vec3 locals_cord(vec3 v, vec3 u1, vec3 u2){
    // return cords of v in local base (u1,u2)
    float a = dot(u2,u2) * dot(v,u1) - dot(u1,u2) * dot(v,u2);
    float b = dot(u1,u1) * dot(v,u2) - dot(u1,u2) * dot(v,u1);
    float k = dot(u1,u1) * dot(u2,u2) - dot(u1,u2) * dot(u2,u1);
    return vec3(a,b,k) ;
}


// Return the reflection of v with respect to the normal n
vec3 reflection(vec3 v, vec3 n){
    return v - 2.0 * dot(v,n) * n / dot(n,n) ;
}

vec3 refraction(vec3 v_, vec3 n, float r){
    vec3 v = normalize(v_);

    float a = dot(n,n) ;
    float b = 2 * r * dot(n,v) ;
    float c = r * r * dot(v,v) - 1 ;

    float d = sqrt(b * b - 4 * a * c);
    if ( dot(v,n) < 0) d*= -1 ;
    float k = (- b + d ) / (2 * a) ;

    return k * n + r * v ;
}

// Return the projection of A on the line L
vec3 projection(vec3 A, Line L){
    // let H be the projection of A on L

    vec3 v = normalize(L.v);
    float d = dot(A-L.point,v) ;
    
    return L.point + d * v;
}

////////////////////////////////////////////////////////////////////////


Line get_intersection(Line L, Sphere S){
    // Return the intersection and the normal of the sphere
    float d = dot(S.center - L.point,normalize(L.v)) ;
    if (d <= 0.0) return Line(vec3(0.),vec3(0.));
    vec3 H = L.point + d * normalize(L.v);
    vec3 u = H - S.center;
    
    float l = S.radius * S.radius - dot(u,u);

    // H isnt in the sphere
    if (l < 0.) return Line(H,vec3(0.));
    
    // H is in the sphere
    // We calculate the intersection based on H
    vec3 I = H - sqrt(l) * normalize(L.v);
    vec3 n = I - S.center;

    return Line(I, n) ;
}

Line get_intersection(Line L, Cylinder C){

    vec3 v1 = cross(L.v,C.v);
    vec3 v2 = cross(L.point - C.origin, C.v);

    
    float a = norm2(v1) ;
    float b = dot( v1, v2 );
    float c = norm2( v2 ) - r * r * norm2(C.v) ;

    float d = b*b - a*c ;

    if (d < 0.) return Line(vec3(0.),vec3(0.));

    float lambda = - (b + sqrt(d)) / a ;
    if (lambda < 0.0) lambda = - (b - sqrt(d)) / a ;
    
    return L.point - lambda * L.v ;
}

Line get_intersection(Line L, Plane T, int type){

    vec3 n = cross(T.u1,T.u2);
    
    
    if (dot(n,L.v) == 0.) return Line(vec3(0.),vec3(0.));

    // The line intersect the plane of the plane
    float lambda = dot(n, T.origin - L.point) / dot(n,L.v);
    if (lambda <= 0.0) return Line(vec3(0.),vec3(0.));
    // H is the intersection
    vec3 H = L.point + lambda * L.v ;
    vec3 v = H - T.origin ;

    vec3 normal = dot(L.v,n) > 0.0 ? -n : n ;

    if (type == 1) return Line(H, normal);

    // a and b are coordinate relative to the local base of the plane (T.u1, T.u2)
    vec3 coordinate = locals_cord(v,T.u1,T.u2);
    if (coordinate.z == 0.0) return Line(vec3(0.),vec3(0.));;
    float a = coordinate.x / coordinate.z;
    float b = coordinate.y / coordinate.z; 
    

    // The intersection is in the rectangle
    if (type == 2){
        if (0. <= a && a <= 1. && 0. <= b && b <= 1.){
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
    if (O.type == 0) return get_intersection(L, O.sphere);
    return get_intersection(L, O.plane, O.type);
}

Object get_intersection(Line Ray){
    float min_len = -1. ;
    Object best_obj ;

    for (int j = 0 ; j < nb_object ; j++){
        Line Intersection = get_intersection(Ray,objects[j]);
        float dist = dot(Intersection.point - Ray.point,Intersection.point - Ray.point);
        if (
            Intersection.v != vec3(0.) && 
            dist != 0 &&
            ( min_len == -1. || dist < min_len )
        ){
            min_len = dist ;
            best_obj = objects[j] ;
        }
    }
    return min_len == -1 ? null_Object() : best_obj ;
}


////////////////////////////////////////////////////////////////////////


vec4 calc_color(Object Obj, Line Normal, float dist_from_origin){

    vec3 col = vec3(0.);
    vec3 obj_col ;
    float n_col = 0.0;

    if ( Obj.color == vec3(0.0, 0.0, 0.01)){ // Color with image
        
        vec3 v = Normal.point - Obj.plane.origin ;
        vec3 pos = locals_cord(v,Obj.plane.u1,Obj.plane.u2);
        obj_col = texture2D(u_tex0,vec2(pos.x / pos.z, 1-pos.y / pos.z) ).rgb;
        // * u_tex0Resolution.xy

    } else if ( Obj.color == vec3(.5,.5,.5)) { // Color the floor
        vec3 H = Normal.point ;
        bool cond = fract(H.x)-fract(H.z) < 0.0;
        obj_col = cond ? vec3(0.25) : vec3(0.75) ;
    } else { // base color
        obj_col = Obj.color ;
    }


    if (Obj.is_light){
        return vec4(Obj.color,-1.0);
    } else if (light_type != 0) { // Add light
        Object Light = objects[0];

        vec3 center = Light.type == 0 ? Light.sphere.center : Light.plane.origin ;
        vec3 light_dir = center - Normal.point ;
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
            float dist = (sqrt(dot(light_dir,light_dir)));
            float coef ;
            coef = 4 * angle / dist ;
            // float coef = 10 / ( dist_from_origin * dist_from_origin) ;

            // col += (1+angle)/2 * (is_illuminate ? vec3(1.0,0,0) :  vec3(0,0,1.0));
            // n_col += (1+angle)/2 ;
            // col += coef * (is_illuminate ? Light.color : vec3(0.0));
            // n_col += coef ;

            if (light_type == 2) col += coef * Light.color * obj_col;
            else col += coef * Light.color;
            n_col += coef ;
        } else {
            col += vec3(0);
            n_col += 1 ;
        }
        
    }

    if (light_type != 2) { // Ambient light
        col += ( 1 - Obj.mirror ) * obj_col ;
        n_col += ( 1 - Obj.mirror ) ;
        
    }
    

    // if (Obj.mirror > 0.){
    //     // col /= objects[j].mirror ;
    //     // n_col /= objects[j].mirror ;
    //     // col *= .5 ;
    //     // n_col *= .5 ;
    //     // Reflection
    //     // Ray = Line(N.point, reflection(N.v,Ray.v));
    // } 
    // else {
    //     Diffusion
    //     vec3 random_dir = vec3(rand(N.point.xy),rand(N.point.xz),rand(N.point.yz));
    //     if (dot(random_dir,N.v) < 0.) random_dir = -random_dir;
    //     Ray = Line(N.point, random_dir);
    // }
    //break;
    return vec4(col,n_col);
}

vec3 draw(Line Ray, Sphere Light){
    vec3 col = vec3(0.);
    float n_col = 0;
    float col_weight = 1;
    float dist_from_origin = 0;


    const int n_reflection = 3;

    for (int i = 0; i < n_reflection; i++){

        Object best_obj = get_intersection(Ray);

        if (best_obj.type == -1) break;

        Line N = get_intersection(Ray,best_obj);
        dist_from_origin += distance(N.point, Ray.point) ;

        vec4 res = calc_color(best_obj, N, dist_from_origin );
        
        if ( res.w == -1.0) {
            col += col_weight * res.xyz;
            n_col += col_weight;
            break ;
        }

        col += col_weight * res.w * res.xyz ;
        n_col += col_weight * res.w ;
        col_weight *= best_obj.mirror ;
        
        Ray.point = N.point ;
        if (best_obj.refrac != 1.0 ){
            Ray.v = refraction(Ray.v,N.v,best_obj.refrac) ;
        } else {
            Ray.v = reflection(Ray.v,N.v);
        }
    }


    if (n_col != 0.){
        return col / n_col ;
    } else if (n_col == 0. && light_type == 2){
        return vec3(0.0);
    } else {
        return vec3(0.0, 0.0, 0.3);
    }
}


void main() {
    float size = max(u_resolution.x,u_resolution.y);
    vec2 st =  1.0 * (gl_FragCoord.xy / vec2(size,size) - .5);
    vec2 mouse = 3 * (u_mouse / u_resolution - .5) ;
	
    

    float t1 = st.x - mouse.x ;
    float t2 = st.y + mouse.y ;
    mat3 M = rot(0,-2*mouse.x,mouse.y);
    // vec3 B = vec3(sin(t1)*cos(t2),sin(t2),cos(t1));
    // vec3 A = vec3(position.x,0,position.y-3) ;
    vec3 B = vec3(st.x,st.y,1.0);
    vec3 A = vec3(position.x,0,position.y-3.0);
    Line Ray = Line(
        A,
        M * B
    );
    
    Sphere Mouse3d = Sphere(
        // 6.0 *  (vec3( mouse.x , -mouse.y ,sin(u_time * 0.0))), .5
        6.0 *  (vec3( 1.0 , -1.0 ,sin(u_time * 0.0))), .5
    );

    vec3 col = draw(Ray,Mouse3d);

    gl_FragColor = vec4(col, 1.0);
   
}


