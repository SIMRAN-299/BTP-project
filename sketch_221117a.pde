  
float h=0;
float m=2;
float x_shift=300;
float y_shift=250;
/*
  float m=22;
float x_shift=140;
float y_shift=-20;// american points
  
  float h=0;
float m=1;
float x_shift=500;
float y_shift=600;//for germany railway

float m=2;
float x_shift=300;
float y_shift=250;//berlin shop


float h=0;
float m=12;
float x_shift=200;
float y_shift=0;//american cities with abrivated names

*/
String[] lines; 
int check=0;

void draw() {
  
  if(check==1){
     lines = loadStrings("output/Berlin_shop_output.txt");
 
  println("there are " + lines.length + " lines");
  for (int i = 0; i < lines.length; i++) {
    String numbers = lines[i];
    if(i==0){
      float[] nums = float(split(numbers, ' '));
      h=nums[0];
    }
    else{
    float[] nums = float(split(numbers, ' '));
    float x=nums[0]+x_shift;
    float y=nums[1]+y_shift;
    float l=nums[2];
    fill(0,127,0);
    rect(x*m,y*m,l*m,h*m);
    }
    
    
    
  }
  }
}
void mouseClicked() {
  check=check+1;
  if(check==2){
    check=0;
  }
}


void setup() {
  
  size(1800, 950);
  background(255);
  lines = loadStrings("testCases/Berlin_shop.txt");
  println("there are " + lines.length + " lines");
  for (int i = 0; i < lines.length; i++) {
    String numbers = lines[i];
   
    float[] nums = float(split(numbers, ' '));
    float x=nums[0]+x_shift;
    float y=nums[1]+y_shift;
    //float l=nums[2];
    fill(0,0,0);
    ellipse(x*m,y*m,2,2);
    //ellipse(x*m,y*m,l*m,h*m);
    
  }
  
 
  
  
}

