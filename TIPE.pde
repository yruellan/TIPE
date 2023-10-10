PShader shader;
String file ;

boolean sky_is_day ;

float x,y ;
float light_h ;
int light_type ;
int picture_n ;
boolean draw_type ;
boolean take_picture ;

void setup() {
  size(800, 800, P2D);
  noStroke();
  println("run Test.pde");
  
  
  file = "shader2.glsl" ;
  shader = loadShader(file);
  
  frameRate(30);
  
  x = -5.0 ;
  y = 0.0 ;
  light_h = 4 ;
  light_type = 0 ;
  picture_n = 0 ;
  draw_type = false ;
  take_picture = false ;
  sky_is_day = true ;
  
  setup_table();
  redraw();
  noLoop();
}


void draw() {
  
  shader.set("u_resolution", float(width), float(height));
  shader.set("u_mouse", float(mouseX), float(mouseY));
  shader.set("position", x,y);
  shader.set("light_height", light_h);
  shader.set("light_type", light_type);
  shader.set("draw_type",draw_type); 
  
  //shader.set("u_time", millis() / 1000.0);
  //shader.set("u_tex0Resolution",1024,1024);
  
  shader.set("u_tex0",loadImage("content/hokusai.jpg")); 
  shader.set("earth",loadImage("content/earth.jpg")); 
  shader.set("sun",loadImage("content/sun.jpg")); 
  shader.set("ground",loadImage("content/ground.jpg")); 
  if (sky_is_day) shader.set("sky",loadImage("content/sky.jpg")); 
  else shader.set("sky",loadImage("content/night_sky.jpg")); 
  
  shader(shader);
  rect(0,0,width,height);
  
  if (take_picture){
    save("pictures/picture"+picture_n+".jpeg");
    picture_n++ ;
    println("picture", picture_n);
    take_picture = false ;
  }
  
}

void keyPressed(){
  float ds = 1 ;
  float dh = 1 ;
  //println("key",keyCode,key);

  if ( keyCode >= 37 && keyCode <= 40){
    PVector spe = new PVector(0,0);
    if      (keyCode == 37){ spe = PVector.fromAngle(PI);}
    else if (keyCode == 38){ spe = PVector.fromAngle(HALF_PI);}
    else if (keyCode == 39){ spe = PVector.fromAngle(0);}
    else if (keyCode == 40){ spe = PVector.fromAngle(-HALF_PI);}
    float t = float(mouseX)/width - .5 ;
    //println("rot",t,2*t,TWO_PI*t);
    spe.rotate(-TWO_PI*t);
    x += spe.x * ds ; y += spe.y * ds ;
  }
  else if (keyCode == 46){ light_h-= dh ;} // :
  else if (keyCode == 47){ light_h+= dh ;} // =
  else if (keyCode == 77){ // ,
    draw_type = !draw_type;
    println("change draw",draw_type);
  } 
  else if (keyCode == 44){
    light_type = (light_type+1)%3 ; // ;
  }
  else if (keyCode == 78){ // n
     sky_is_day= !sky_is_day;
  }
  else if (keyCode == 10){
    take_picture = true ;
    //save("pictures/picture"+picture_n+".jpeg");
    //picture_n++ ;
    //println("picture", picture_n);
  }
  else if (keyCode == 27){
    // escape
    exit();
  }
  else{
    print("Reaload Shader  ",key,": ",keyCode);
    shader = loadShader(file);
    println("  ;");
  }
  redraw();
}

void mouseMoved(){
  redraw();
}
