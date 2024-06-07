
// https://www.youtube.com/watch?v=BFld4EBO2RE Video Painting a Landscape with Maths (IQ)

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


//////////////////////////////////////////////////////////////////////

#define PI 3.1415926535897932384626433832795
#define TWO_PI 6.283185307179586476925286766559


//////////////////////////////////////////////////////////////////////

struct Line {
    vec3 origin ;
    vec3 v ;
};

//////////////////////////////////////////////////////////////////////

mat2 rot(float t){
    return mat2(
        cos(t),sin(t),
        -sin(t),cos(t)
    );
}

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

float random(float x){
    return fract(sin(x)*10000.0) ;
}

//////////////////////////////////////////////////////////////////////

float noise(vec2 p){
    return sin(p.x) * sin(p.y);
}

float S(float s){
    return 3*s*s - 2*s*s*s ;
}

float A(float i, float j){
    // return 1 ;
    // vec2 p = fract(vec2(i,j*2.1827)/PI) ;
    // return 2*fract( p.x*p.y*(p.x+p.y) )-1 ;
    vec2 p = 50.0*fract( vec2(i,j)*0.3183099 );
    return fract( p.x*p.y*(p.x+p.y) ); 
}

float N(vec2 p){
    float scale = 1 ;
    float i = floor(p.x);
    float j = floor(p.y);

    float a = A(i,j) ;
    float b = A(i+1,j) ;
    float c = A(i,j+1) ;
    float d = A(i+1,j+1) ;

    return a + (b-a)*S(p.x - i) + (c-a)*S(p.y-j) + (a-b-c+d)*S(p.x-i)*S(p.y-j) ;
}

vec4 octaves(vec2 p){
    float h = 0 ;
    float k = 1 ;

    vec3 normal = vec3(0,0,0) ;
    normal.y = 1;

    mat2 R = mat2( .8, .6, -.6, .8) ;
    float a = 1.7 ;

    for (int i = 0 ; i < 8 ; i++){
        h += N(p) / k;

        normal.x += (N(p)-N(p+vec2(0.01,0))) / (k* 0.01);
        normal.z += (N(p)-N(p+vec2(0,0.01))) / (k* 0.01);

        // p += 5*vec2(random(i), random(i+34.227));
        p = a * R * p ;
        k *= a;
    }
    return vec4(normalize(normal),h) ;
}

float shadow(vec3 pos, vec3 sun_dir){

    float R_min = 1.0 ;

    float ds = 0.01 ;
    vec3 p = pos ;;

    for (int i = 0 ; i < 500 ; i++){
        p += ds * sun_dir ;
        float h = octaves(p.xz).w;

        // R_min = min(R_min,32.0 * (p.y - h)/(float(i) * 0.1)) ;
        if (p.y< h) {
            R_min = 0 ;
            break;
        }
    }

    return clamp(0,1,R_min) ;
}

//////////////////////////////////////////////////////////////////////

vec4 raymarch(vec3 p, vec3 dir){
    const int N = 700 ;
    float ds = 0.01 ;

    for ( int i = 0 ; i<N ; i++){
        // p += ds * dir ;
        p += ds * dir * float(i) ;
        vec4 oct = octaves(p.xz/3);
        float h = oct.w;

        if (p.y < h) {
            // normal = ;
            return vec4(oct.xyz,float(i)/N);
        }
    }

    return vec4(vec3(0),1.0) ;
}

vec3 draw(Line Ray){

    const int N = 700 ;
    float ds = .01 ;

    float sun_angle = variable*TWO_PI/10;
    vec3 sun_dir = normalize(vec3(cos(sun_angle),1,sin(sun_angle)));
    
    vec4 Obj = raymarch(Ray.origin, Ray.v) ;
    float r = Obj.w ;
    vec3 normal = Obj.xyz;

    if (r == 1.0){ // sky
        vec3 col = vec3(0.0, 0.749, 1.0) - Ray.v.y * 0.5 ;

        float k = clamp(0,1,dot(sun_dir, Ray.v));
        k = k*k * .4 ;
        col = vec3(1.0, 1.0, 0.0) * k + (1-k) * col;

        // Cloud : 13:56
        // float lambda = raymarch(0.1*Ray.origin+vec3(0,20,0),Ray.v*vec3(1,-1,1)).w;
        // lambda = sqrt(lambda) ;
        // col = col * (1-lambda) + lambda * vec3(1.0, 1.0, 1.0) ;

        return col ;
    }

    vec3 col = vec3(r);
    col.r = clamp(0,1, dot(sun_dir, normal)) * (1-r)/2;
    col.g = clamp(0,1, dot(sun_dir, normal)) * (1-r)/2;
    // col = (.5+.5*shadow(vec3(p.x,h,p.y), sun_dir)) * col ;

    float k = 0;
    // k = pow(2,-r*3);
    k = 1-r ;

    // return col ;
    return col * k + (1-k) * vec3(0.5,0.5,0.5);
}

void main() {
    float size = max(u_resolution.x,u_resolution.y);
    vec2 st =  1.0 * (gl_FragCoord.xy / vec2(size,size) - .5);

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