PShader shader;
String file ;



float x,y ;
float angleX, angleY ;

float light_h ;
int light_type ;
boolean sky_is_day ;
int draw_type ;

boolean scene1 = true ;
boolean scene2 = false ;
boolean scene3 = false ;
boolean bases_vectors = false ; 

int picture_n ;
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
  draw_type = 0 ;
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
    print("Reload");
    flush();
    shader.set("u_mouse", width/2.0+angleX, height/2.0+angleY);
    shader.set("light_height", light_h);
    shader.set("draw_type",draw_type); 
    shader.set("light_type", light_type);
    shader.set("position", x,y);
    
    shader.set("scene1",scene1);
    shader.set("scene2",scene2);
    shader.set("scene3",scene3);
    shader.set("bases_vectors",bases_vectors);
    
    shader.set("u_resolution", float(width), float(height));
    shader.set("hokusai",loadImage("content/hokusai.jpg")); 
    shader.set("earth",loadImage("content/earth.jpg")); 
    shader.set("sun",loadImage("content/sun.jpg")); 
    shader.set("ground",loadImage("content/ground.jpg")); 
    shader.set("galaxy",loadImage("content/galaxy.png")); 
    shader.set("moon",loadImage("content/moon.jpg")); 
    shader.set("sagittarius_A",loadImage("content/sagittarius_A.jpg")); 
    
    if (sky_is_day) shader.set("sky",loadImage("content/sky.jpg")); 
    else shader.set("sky",loadImage("content/night_sky.jpg"));
    shader(shader);
    reload_shader = false ;
    println("\tShader Reloaded");
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
    //float t = (mouseX+angleX)/width - .5 ;
    float t = angleX/width ;
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
    
    shader.set("u_mouse", width/2.0+angleX, height/2.0+angleY);
  }
  
  else if (keyCode == 49){ // 1&
    scene1 = !scene1 ;
    shader.set("scene1",scene1);
  }
  else if (keyCode == 50){ // 2é
    scene2 = !scene2 ;
    shader.set("scene2",scene2);
  }
  else if (keyCode == 51){ // 3"
    scene3 = !scene3 ;
    shader.set("scene3",scene3);
  }
  else if (keyCode == 52){ // 4'
    bases_vectors = !bases_vectors ;
    shader.set("bases_vectors",bases_vectors);
  }
  
  else if (keyCode == 46){  // :
    light_h-= dh ; shader.set("light_height", light_h);}
  else if (keyCode == 47){  // =
    light_h+= dh ; shader.set("light_height", light_h);}
  else if (keyCode == 77){  // ,
    draw_type = (draw_type+1)%3;
    println("change draw type : ",draw_type);
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
  else if (keyCode == 10){ // ENTER ?
    take_picture = true ;
    //save("pictures/picture"+picture_n+".jpeg");
    //picture_n++ ;
    //println("picture", picture_n);
  }
  else if (keyCode == 27){ // ESC
    // escape
    exit();
  }
  else if (keyCode == 32){ // SPACE
    shader = loadShader(file);
    reload_shader = true ;
  } else {
    println("Key : ",key,": ",keyCode);
  }
  redraw();
}

void mouseMoved(){
  //shader.set("u_mouse", mouseX+angleX, mouseY+angleY);
  //redraw();
}
