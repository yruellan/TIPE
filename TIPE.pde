PShader shader;
String file ;

boolean sky_is_day ;

float x,y ;
float angleX, angleY ;
float light_h ;
int light_type ;
int picture_n ;
boolean draw_type ;
boolean take_picture ;
boolean reload_shader ;

void setup() {
  size(800, 800, P2D);
  //size(600, 600, P2D);
  noStroke();
  println("run Test.pde");
  
  file = "shader.glsl" ;
  shader = loadShader(file);
  
  frameRate(30);
  
  x = 0.0 ;
  y = -5.0 ;
  angleX = 0 ;
  angleY = 0 ;
  light_h = 4 ;
  light_type = 0 ;
  picture_n = 0 ;
  draw_type = false ;
  take_picture = false ;
  reload_shader = true ;
  sky_is_day = true ;
  
  setup_table();
  
  noLoop();
  redraw();
 
 //https://forum.processing.org/one/topic/run-code-on-exit.html
 
}


void draw() {
  
  rect(0,0,width,height);
  
  if (take_picture){
    //save("pictures/picture"+picture_n+".jpeg"); // Bad quality
    save("pictures/picture"+picture_n+".png");
    picture_n++ ;
    println("picture", picture_n);
    take_picture = false ;
  }
  if (reload_shader){
    shader.set("u_mouse", float(width/2), float(height/2));
    shader.set("light_height", light_h);
    shader.set("draw_type",draw_type); 
    shader.set("light_type", light_type);
    shader.set("position", x,y);
    
    shader.set("u_resolution", float(width), float(height));
    shader.set("hokusai",loadImage("content/hokusai.jpg")); 
    shader.set("earth",loadImage("content/earth.jpg")); 
    shader.set("sun",loadImage("content/sun.jpg")); 
    shader.set("ground",loadImage("content/ground.jpg")); 
    shader.set("galaxy",loadImage("content/galaxy.png")); 
    
    if (sky_is_day) shader.set("sky",loadImage("content/sky.jpg")); 
    else shader.set("sky",loadImage("content/night_sky.jpg"));
    shader(shader);
    reload_shader = false ;
    println("Shader Reloaded");
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
    float t = (mouseX+angleX)/width - .5 ;
    //println("rot",t,2*t,TWO_PI*t);
    spe.rotate(-TWO_PI*t);
    x += spe.x * ds ; y += spe.y * ds ;
    shader.set("position", x,y);
  }
  else if (key == 'z' || key == 'q' ||key == 's' ||key == 'd'){
    if      (key == 'z'){ angleY -= 6;}
    else if (key == 'q'){ angleX -= 6;}
    else if (key == 's'){ angleY += 6;}
    else if (key == 'd'){ angleX += 6;}
    
    shader.set("u_mouse", mouseX+angleX, mouseY+angleY);
  }
  
  else if (keyCode == 46){  // :
    light_h-= dh ; shader.set("light_height", light_h);}
  else if (keyCode == 47){  // =
    light_h+= dh ; shader.set("light_height", light_h);}
  else if (keyCode == 77){  // ,
    draw_type = !draw_type;
    println("change draw",draw_type);
    shader.set("draw_type",draw_type); 
  } 
  else if (keyCode == 44){  // ;
    light_type = (light_type+1)%3 ; 
    shader.set("light_type", light_type);
  }
  else if (keyCode == 78){ // n
     sky_is_day= !sky_is_day;
     if (sky_is_day) shader.set("sky",loadImage("content/sky.jpg")); 
     else shader.set("sky",loadImage("content/night_sky.jpg")); 
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
    println("Reload Shader  ",key,": ",keyCode);
    shader = loadShader(file);
    reload_shader = true ;
  }
  redraw();
}

void mouseMoved(){
  shader.set("u_mouse", mouseX+angleX, mouseY+angleY);
  redraw();
}
