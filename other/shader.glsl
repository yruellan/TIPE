// Author:
// Title:
// https://raytracing.github.io/
// https://github.com/RayTracing/raytracing.github.io/
// https://images.math.cnrs.fr/Un-arbre-pythagoricien.html
// https://www.josleys.com/
// https://cal.cs.umbc.edu/Courses/CMSC435-F15/Slides/raytrace.pdf
// https://www.cs.cornell.edu/courses/cs4620/2014fa/lectures/05rt-shading.pdf

// http://www.joelsornette.fr/Archives/exotypes/exotype26.pdf

// #ifdef GL_ES
// precision mediump float;
precision highp float;
// #endif

uniform sampler2D u_tex0;
uniform vec2 u_tex0Resolution;

uniform vec2 u_resolution;
uniform vec2 u_mouse;
uniform float u_time;
uniform vec2 position;
uniform float light_height;
uniform int light_type; // 0: no light, 1: light+ambiente, 2: only light
uniform bool draw_type ;


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
    // -1 error, 0 sphere, 1 cylinder, 2 plane, 3 rect, 4 circle, 5+ triangle
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
        -1,                                         // type
        vec3(0.0),                                // color
        0.0,                                        // mirror
        1.0,                                        // refrac
        false                                       // is_light
    ) ;
}

Object new_Object(Sphere S, vec3 color, float mirror){
    Object O = null_Object();
    O.sphere = S;
    O.type = 0;
    O.color = color;
    O.mirror = mirror;
    return O;
}
Object new_Object(Cylinder C, vec3 color, float mirror){
    Object O = null_Object();
    O.cylinder = C;
    O.type = 1;
    O.color = color;
    O.mirror = mirror;
    return O;
}
Object new_Object(Plane P, vec3 color, float mirror){
    Object O = null_Object();
    O.plane = P;
    O.type = 2;
    O.color = color;
    O.mirror = mirror;
    return O;
}

Object new_Rectangle(Plane P, vec3 color, float mirror){
    Object O = new_Object(P, color, mirror);
    O.type = 3;
    return O;
}
Object new_Plane(Plane P, vec3 color, float mirror){
    Object O = new_Object(P, color, mirror);
    O.type = 2;
    return O;
}
Object new_typed_Plan(Plane P, vec3 color, float mirror, int type){
    Object O = new_Object(P, color, mirror);
    O.type = type;
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

#define nb_object 7


Object objects[nb_object] = Object[nb_object](
    Light(new_Object(Sphere(vec3(0,light_height,0),.2),vec3(0.9, 0.85, 1.0),0.0)),
    //Light(new_Rectangle(
    //    Plane(vec3(-1,light_height,1),vec3(2,0,0),vec3(0,0,2)),
    //    vec3(0.9, 0.85, 1.0),0
    //)),
    //new_Object(Sphere(vec3(1,0,-0.5),.5),vec3(0.0, 0.7, 0.0),0.0),
    new_Object(Sphere(vec3(2,0.1,5),2),vec3(0.8824, 1.0, 0.0),0.8),
    new_Object(Sphere(vec3(2,0.1,-5),2),vec3(1.0, 0.0, 0.0),0.8),
    //new_Object(Sphere(vec3(-1,0,-0.5),.5),vec3(0.99, 0.0, 0.0),0.0),
    
    // Floor :
    new_Object( 
        Plane(vec3(0,-2,0),vec3(2,0,0),vec3(0,0,2)),
        vec3(0.5, 0.5, 0.5),0.0
    ),

    // Mirrors :

    new_typed_Plan(
       Plane(vec3(-10,-1,10),vec3(20,0,0),vec3(0,6,0)),
       vec3(0.4, 1.0, 0.0),0.0,2
    ),
    // new_Object(
    //     Plane(vec3(15,0,0),vec3(0,0,2),vec3(0,2,0)),
    //     vec3(0.0039, 0.0039, 0.0078),0.8
    // ),

    // Objects

    // new_typed_Plan(
    //     Plane(vec3(5,0.1,5),vec3(0,2,0),vec3(2,0,2)),
    //     vec3(0,1,0), 0.0, 4
    // ),
    // new_Rectangle(
    //     Plane(vec3(3,-.5,-5),vec3(0,4,0),vec3(-8,0,-2)),
    //     vec3(1,0,0), .9
    // ),
    new_Object(Cylinder(vec3(-8,-2,5),vec3(0,5,1),1),vec3(0.0, 0.7, 1.0),0.),
    // new_Object(Cylinder(vec3(-6,-2,-1),vec3(0,6,6),.3),vec3(0.8, 1.0, 0.0),0),

    // picture :

    new_Rectangle(
        Plane(vec3(4,1,0),vec3(-1,0,-3),vec3(0,2,0)),
        vec3(0,0,0.01), 0.0
    )

    // Objects for refraction :

    // new_Object(Sphere(vec3(-10,-1,5),2),vec3(1,0,1),1)
    // Refractor(new_Object(Sphere(vec3(-10,-1,10),2),vec3(1,0,1),1),1.1)
    // Refractor(new_Object(Cylinder(vec3(-6,-2,0),vec3(0,2.5,0),3),vec3(1, 0, 0),1), 1.4)
    // Refractor(new_Rectangle(Plane(vec3(-5,-1,0),vec3(0,3,0),vec3(0,0,3)),vec3(1,0,1),.8),1.4),
    // Refractor(new_Rectangle(Plane(vec3(-5,-1,0),vec3(0,3,0),vec3(-3,0,0)),vec3(1,0,1),.8),1.4)
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


// Return the reflexion of v with respect to the normal n
vec3 reflexion(vec3 v, vec3 n){
    return v - 2.0 * dot(v,n) * n / dot(n,n) ;
}

vec3 refraction(vec3 v_, vec3 n, float r){
    ///////////////////
    //  r = n1 / n2  //
    ///////////////////
    vec3 v = normalize(v_);

    float a = dot(n,n) ;
    float b = r * dot(n,v) ;
    float c = r * r * dot(v,v) - 1 ;

    if (b * b - a * c < 0) return vec3(0.);
    // if (b * b - a * c < 0) return reflexion(v,n);
    float d = sqrt(b * b - a * c);
    if ( dot(v,n) < 0) d*= -1 ;
    float k = (- b + d ) / a ;

    return k * n + r * v ;
}

vec2 taux_ref(float i1, float i2,float r){
    // https://claude-gimenes.fr/physique/propagation-des-ondes-electro-magnetiques/-ii-reflexion-et-refraction-en-milieux-isotropes-1
    // https://fr.wikipedia.org/wiki/Coefficient_de_Fresnel
    // http://www.joelsornette.fr/Archives/exotypes/exotype26.pdf

    // float r = (i1 - i2) / (i1 + i2) ;
    // float t = 2 * i1 / (i1 + i2) ;

    // float r = tan(i2-i1) / tan(i1+i2) ;
    // float t = 2 * sin(i2) * cos(i1) / (sin(i1+i2) * cos(i2-i1)) ;

    // float r = (cos(i1) - r*cos(i2)) / (cos(i1) + r*cos(i2));
    // float t = 2 * cos(i1) / (cos(i1) + r*cos(i2));

    if (r < 1 && sqrt(1-r*r) >= 1) return vec2(1,0);

    float n1 = 1 ;
    float n2 = r * n1 ;

    float R = (n2 * cos(i1) - n1 * cos(i2))**2 / (n2 * cos(i1) + n1 * cos(i2))**2 ;
    float T = 4*n1*n2*cos(i1)*cos(i2)) / (n2 * cos(i1) + n1 * cos(i2))**2 ;

    return vec2(R,T);
}

vec2 taux_ref_2(vec3 u, vec3 v_t, vec3 n, float r){
    // http://www.joelsornette.fr/Archives/exotypes/exotype26.pdf
    // u Light vector, v_t transmitted vector, n normal vector, r quotient of refraction

    if (r < 1 && sqrt(1-r*r) >= 1) return vec2(1,0);

    float n1 = 1 ;
    float n2 = r * n1 ;

    float R = (n2 * dot(u,n) - n1 * dot(v_t,n))**2 / (n2 * dot(u,n) + n1 * dot(v_t,n))**2 ;
    float T = 4*n1*n2*dot(u,n)*dot(v_t,n)) / (n2 * dot(u,n) + n1 * dot(v_t,n))**2 ;

    return vec2(R,T);
}

////////////////////////////////////////////////////////////////////////


Line get_intersection(Line L, Sphere S){
    // Return the intersection and the normal of the sphere
    float d = dot(S.center - L.origin,normalize(L.v)) ;
    if (d <= 0.0) return Line(vec3(0.),vec3(0.));
    vec3 H = L.origin + d * normalize(L.v);
    vec3 u = H - S.center;
    
    float l = S.radius * S.radius - dot(u,u);

    // H isnt in the sphere
    if (l < 0.) return Line(H,vec3(0.));
    
    // H is in the sphere
    // We calculate the intersection based on H
    
    vec3 I = H ;
    if (d - sqrt(l) < 0) I += sqrt(l) * normalize(L.v);
    else I -= sqrt(l) * normalize(L.v);
    // if (dot(I-L.origin, L.v) < 0) I = H - sqrt(l) * normalize(L.v);

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
    
    
    if (dot(n,L.v) == 0.) return Line(vec3(0.),vec3(0.));

    // The line intersect the plane of the plane
    float lambda = dot(n, T.origin - L.origin) / dot(n,L.v);
    if (lambda <= 0.0) return Line(vec3(0.),vec3(0.));
    // H is the intersection
    vec3 H = L.origin + lambda * L.v ;
    vec3 v = H - T.origin ;

    vec3 normal = dot(L.v,n) > 0.0 ? -n : n ;

    if (type == 2) return Line(H, normal);

    // a and b are coordinate relative to the local base of the plane (T.u1, T.u2)
    vec3 coordinate = locals_cord(v,T.u1,T.u2);
    if (coordinate.z == 0.0) return Line(vec3(0.),vec3(0.));;
    float a = coordinate.x / coordinate.z;
    float b = coordinate.y / coordinate.z; 
    

    // The intersection is in the rectangle
    if (type == 3){
        if (0. <= a && a <= 1. && 0. <= b && b <= 1.){
            return Line(H, normal);
        } else {
            return Line(vec3(0.),vec3(0.));
        }
    }

    // The intersection is in the circle
    if (type == 4){
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
    if (O.type == 0) return get_intersection(L, O.sphere);
    if (O.type == 1) return get_intersection(L, O.cylinder);
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
        
        vec3 v = Normal.origin - Obj.plane.origin ;
        vec3 pos = locals_cord(v,Obj.plane.u1,Obj.plane.u2);
        obj_col = texture2D(u_tex0,vec2(pos.x / pos.z, 1-pos.y / pos.z) ).rgb;
        // * u_tex0Resolution.xy

    } else if ( Obj.color == vec3(.5,.5,.5)) { // Color the floor
        vec3 H = Normal.origin ;
        bool cond = fract(H.x)-fract(H.z) < 0.0;
        obj_col = cond ? vec3(0.25) : vec3(0.75) ;
    } else { // base color
        obj_col = Obj.color ;
        // obj_col = normalize(Normal.v) ;
    }


    if (Obj.is_light){
        return vec4(Obj.color,-1.0);
    } else if (light_type != 0 && Obj.mirror != 1.0) { // Add light
        

        for (int i=0 ; i < nb_object ; i++) {

            Object Light = objects[i];
            if ( ! Light.is_light ){ continue; }

            vec3 center = Light.type == 0 ? Light.sphere.center : Light.plane.origin ;
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
                float dist =norm2(light_dir);
                float coef ;
                coef = angle * 6 / (1.0 + dist) ;
                // coef = 4 * angle*angle / dist ;
                // coef = 4  / dist ;
                // float coef = 10 / ( dist_from_origin * dist_from_origin) ;

                // col += (1+angle)/2 * (is_illuminate ? vec3(1.0,0,0) :  vec3(0,0,1.0));
                // n_col += (1+angle)/2 ;
                // col += coef * (is_illuminate ? Light.color : vec3(0.0));
                // n_col += coef ;
                if (angle != 0){
                    if (light_type == 2) col += coef * Light.color * obj_col;
                    else col += coef * Light.color;
                    n_col += coef ;
                }
            } else {
                col += vec3(0);
                n_col += 1 ;
            }

        }
    }

    if (light_type != 2) { // Ambient light
        col += ( 1 - Obj.mirror ) * obj_col ;
        n_col += ( 1 - Obj.mirror ) ;
        
    }
    
    return vec4(col,n_col);
}

vec3 draw(Line Ray, const int n_reflection){
    vec3 col = vec3(0.);
    float n_col = 0;
    float col_weight = 1;
    float dist_from_origin = 0;


    // const int n_reflection = 4;

    for (int i = 0; i < n_reflection; i++){

        Object best_obj = get_intersection(Ray);

        if (best_obj.type == -1){
            if (light_type == 2) col += col_weight * vec3(0.0) ;
            else col += col_weight * vec3(0.0, 0.0, 0.3) ;
            n_col += col_weight ;
            break ;
        }

        Line N = get_intersection(Ray,best_obj);
        dist_from_origin += distance(N.origin, Ray.origin) ;

        vec4 res = calc_color(best_obj, N, dist_from_origin );
        
        if ( res.w == -1.0) {
            col += col_weight * res.xyz;
            n_col += col_weight;
            break ;
        }
 

        col += col_weight * res.w * res.xyz ;
        n_col += col_weight * res.w ;
        col_weight *= best_obj.mirror ;

        if (col_weight == 0) break ;
        
        if (best_obj.refrac != 1.0 ){
            bool is_reverse ;
            if (best_obj.type == 0) is_reverse = dot(N.origin - best_obj.sphere.center,Ray.v) > 0.0 ;
            else if (best_obj.type == 1) is_reverse = dot(N.v,Ray.v) > 0.0 ;
            // else if (best_obj.type == 1) is_reverse = dot(N.origin - best_obj.cylinder.origin,Ray.v) > 0.0 ;
            else is_reverse = dot(cross(best_obj.plane.u1,best_obj.plane.u2),Ray.v) > 0 ;
            // col += 3 * (is_reverse ? vec3(.5,.1,.1) : vec3(.1,.1,.5)) ;
            // n_col += 3 ;
            // break ;
            float r = is_reverse ? 1/best_obj.refrac : best_obj.refrac ;
            Ray.v = refraction(Ray.v,N.v,r) ;
            if (Ray.v == vec3(0.0)) break ;
            // Ray.v = refraction(Ray.v,N.v,best_obj.refrac) ;
        } else {
            Ray.v = reflexion(Ray.v,N.v);
        }
        Ray.origin = N.origin ;

        // iter loop for calc reflexion or refraction
    }

    return col / n_col ;
}

vec3 draw2(Line Ray, const int n_rays_){
    
    const int n_rays = 4;

    Line Rays[n_rays] ;
    vec3 cols[n_rays] ;
    float coeffs[n_rays] ;
    bool is_init[n_rays] ;

    // Init Array
    for (int i = 0; i < n_rays; i++){
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

    float dist_from_origin = 0 ;
    

    for (int i = 0; i < n_rays; i++){

        if (! is_init[i]) continue ;
        Object best_obj = get_intersection(Rays[i]);

        if (best_obj.type == -1){
            if (light_type == 2) cols[i] =  vec3(0.0) ;
            else cols[i] = vec3(0.0, 0.0, 0.3) ;

            continue ;
        }

        Line N = get_intersection(Ray,best_obj);
        dist_from_origin += distance(N.origin, Ray.origin) ;

        vec4 res = calc_color(best_obj, N, dist_from_origin );
        
        // if ( res.w == -1.0) { // obj is light
        //     cols[i] = res.xyz;

        //     continue ;
        // }
        cols[i] = res.rgb ;
        // cols[i] = res.xyz * res.w ;
        // coeffs[i] *= res.w ;

        // best_obj.mirror != 0.0 && 
        
        if (false && n-2< n_rays && best_obj.refrac != 1){

            bool is_reverse ;
            if (best_obj.type == 0) is_reverse = dot(N.origin - best_obj.sphere.center,Ray.v) > 0.0 ;
            else if (best_obj.type == 1) is_reverse = dot(N.v,Ray.v) > 0.0 ;
            // else if (best_obj.type == 1) is_reverse = dot(N.origin - best_obj.cylinder.origin,Ray.v) > 0.0 ;
            else is_reverse = dot(cross(best_obj.plane.u1,best_obj.plane.u2),Ray.v) > 0 ;
            // col += 3 * (is_reverse ? vec3(.5,.1,.1) : vec3(.1,.1,.5)) ;
            // n_col += 3 ;
            // break ;
            float r = is_reverse ? 1/best_obj.refrac : best_obj.refrac ;

            Rays[n].v = reflexion(Rays[i].v,N.v);
            Rays[n].origin = N.origin ;
            float C = taux_ref_2(u,Rays[n].v,N.v,r) ;
            coeffs[n] *= C.x ;
            is_init[n] = true ;
            n++ ;
            Rays[n].v = refraction(Rays[i].v,N.v);
            Rays[n].origin = N.origin ;
            coeffs[n] *= C.y ;
            is_init[n] = true ;
            n++ ;

        } else if (best_obj.mirror != 0.0 && n-1<n_rays) {
            Rays[n].v = reflexion(Rays[i].v,N.v);
            Rays[n].origin = N.origin ;
            is_init[n] = true ;
            n++ ;
        }
    }

    vec3 col = vec3(0) ;
    float sum_coeffs = 0 ;
    int nb_init = 0 ;

    for (int i = 0; i < n_rays; i++){
        if (is_init[i]){
            col += cols[i] ; 
            sum_coeffs += coeffs[i] ;
            nb_init += 1 ;
        }
    }
    // return vec3(2,1,0.1) / 1 ;
    // return col ;
    // return col / sum_coeffs ;
    // if (nb_init == 0) return vec3(0);
    // if (nb_init == 1) return vec3(1,0,0);
    // if (nb_init == 2) return vec3(0,1,0);
    // if (nb_init == 3) return vec3(0,0,1);
    // if (nb_init == 4) return vec3(0,1,1);
    // if (nb_init == 5) return vec3(1,0,1);
    // if (nb_init == 6) return vec3(0.4);
    // if (nb_init == 7) return vec3(1.0, 0.7, 0.0);
    // return vec3(1);
    if (sum_coeffs == 0) return vec3(0.3, 0.6, 1.0) ;
    if (nb_init == 0) return vec3(0.8039, 0.149, 0.7686) ;
    return col / nb_init ;
    // return col / nb_init ;
}

void main() {
    float size = max(u_resolution.x,u_resolution.y);
    vec2 st =  1.0 * (gl_FragCoord.xy / vec2(size,size) - .5);
    const int n_rays = 4 ;

    // vec3 col = vec3(0);
    // const int calc_by_pixel = 1 ;

    // for (int i = 0; i < calc_by_pixel; i++){
    //     for (int j = 0; j < calc_by_pixel; j++){
    //         vec2 mouse = 3 * ((u_mouse + vec2(i,j)/calc_by_pixel) / u_resolution - .5) ;
    //         mat3 M = rot(0,-2*mouse.x,mouse.y);
    //         vec3 B = vec3(st.x,st.y,0.66);
    //         vec3 A = vec3(position.x,0,position.y);

    //         Line Ray = Line(
    //             A,
    //             normalize(M * B)
    //         );
            
    //         col += (
    //             draw_type ? draw(Ray,n_rays) : draw(Ray,n_rays)
    //         ) / (calc_by_pixel*calc_by_pixel);
    //         // if (draw_type) {col += draw(Ray,n_rays) / (calc_by_pixel*calc_by_pixel);}
    //         // else {col += draw2(Ray,n_rays) / (calc_by_pixel*calc_by_pixel);}
    //     }
    // }

    vec2 mouse = 3 * (u_mouse / u_resolution - .5) ;
    mat3 M = rot(0,-2*mouse.x,mouse.y);
    vec3 B = vec3(st.x,st.y,0.66);
    vec3 A = vec3(position.x,0,position.y);

    Line Ray = Line(
        A,
        normalize(M * B)
    );
    vec3 col = draw_type ? draw(Ray,n_rays) : draw2(Ray,n_rays) ;
    
    gl_FragColor = vec4(col, 1.0);
   
}


