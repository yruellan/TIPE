// just the draw function (the same in shader.glsl)

vec3 draw(Line Ray){
    vec3 col = vec3(0.);
    float n_col = 0;
    float col_weight = 1;
    float dist_from_origin = 0;


    const int n_reflection = 4;

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