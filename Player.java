import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

class Player {
    private final static int MARS_WIDTH = 7000;
    private final static int MARS_HEIGHT = 3000;
    private static double[] landscape;

    public static void main(String args[]) {
        Scanner in = new Scanner(System.in);
        int N = in.nextInt(); // the number of points used to draw the surface of Mars.
        landscape = new double[N * 2];
        for (int i = 0; i < N*2; i+=2) {
            landscape[i]     = in.nextInt() / 4.0; // X coordinate of a surface point. (0 to 6999)
            landscape[i + 1] = 750 - in.nextInt() / 4.0; // Y coordinate of a surface point. By linking all the points together in a sequential fashion, you form the surface of Mars.
        }
        for(double i : landscape){
            System.err.print(i + ", ");
        }
        System.err.println();
        
        while (true) {
            GameState gameState = new GameState(in.nextInt(), in.nextInt(), in.nextInt(),
                    in.nextInt(), in.nextInt(), in.nextInt(), in.nextInt());
            
            System.out.println("-20 3");
        }
    }
    static double[] getNearestLandscapeSegment(double x){
        for (int i = 0; i < landscape.length; i+=2){
            if(landscape[i] >= x) 
                return new double[]{landscape[i-2], landscape[i-1],landscape[i], landscape[i+1]};
        }
        return new double[]{0.0, 0.0, 0.0, 0.0};
    }
    
    static class GeneticAlgorithm{
        public void evaluate(){
            
        }
    }
    
    static class Individual {
        ArrayList<Gene> path;
        
    }
    
    static class Gene{
        int angle;
        int power;
        int time;
        public Gene(int angle, int power, int time) {
            this.angle = angle;
            this.power = power;
            this.time  = time;
        }
    }
    
    static class GameState{
        double x;
        double y;
        double hSpeed;
        double vSpeed;
        int fuel;
        int angle;
        int power;
        
        public GameState(double x, double y, double hSpeed, double vSpeed, int fuel, int angle,
                         int power) {
            this.x      = x;
            this.y      = y;
            this.hSpeed = hSpeed;
            this.vSpeed = vSpeed;
            this.fuel   = fuel;
            this.angle  = angle;
            this.power  = power;
        }

        public GameState(GameState gs) {
            x      = gs.x;
            y      = gs.y;
            hSpeed = gs.hSpeed;
            vSpeed = gs.vSpeed;
            fuel   = gs.fuel;
            angle  = gs.angle;
            power  = gs.power;
        }

        public List<GameState> simulate(int desiredAngle, int desiredPower, int time){
            int dPower = desiredPower - power;
            int dAngle = desiredAngle - angle;
            List<GameState> history = new ArrayList<>();
            history.add(new GameState(this));
            while((Math.abs(dPower)>1 || Math.abs(dAngle)>15) && time>0){
                if(dPower!=0){
                    power += dPower>0? 1:-1;
                    dPower+= dPower>0? -1:1;
                }
                if(dAngle!=0){
                    angle += dAngle>0? Math.min(dAngle,15):Math.max(dAngle, -15);
                    dAngle+= dAngle>0? Math.max(dAngle, -15):Math.min(dAngle, 15);
                }
                move(1);
                time--;
                history.add(new GameState(this));
            }
            if(time>0) {
                move(time);
                history.add(new GameState(this));
            }
            return history;
        }
        
        private void move(int time){
            double radians = Math.toRadians(angle);;
            double xAcc = Math.sin(radians) * power * time;
            double yAcc = (Math.cos(radians) * power - 3.711) * time;
            x += hSpeed - xAcc * 0.5;
            y += vSpeed + yAcc * 0.5;
            hSpeed -= xAcc;
            vSpeed += yAcc;
            fuel -= time*power;
            double[] land = getNearestLandscapeSegment(x);
            isPointBelow(x,y,land);
        }
        
        
        private boolean isPointBelow(double pX,double pY,double[] land){
            double vX = land[2]-land[0];
            double vY = land[3]-land[1];
            double vPx = pX - land[0];
            double vPy = pY - land[1];
            double dot = vX * vPy - vY * vPx;
            return dot>0;
        }
        
        
        
        private boolean isIntersect(double[]traj, double[] land){
            double endAlitude = traj[3];
            if(endAlitude > land[1] && endAlitude > land[3]) return false;
            if (endAlitude < land[1] && endAlitude < land[3]) return true;
            
        }

        private boolean ccw(double x1, double y1, double x2,double y2,double x3,double y3){
            return (y3-y1) * (x2-x1) > (y2-y1) * (x3-x1);
        }
                
        private boolean isIntersect(double x1, double y1, double x2, double y2, double x3,double y3, double x4,double y4){
            return ccw (x1, y1, x3, y3, x4, y4) != ccw(x2, y2, x3,y3, x4, y4) && ccw(x1, y1, x2,y2,x3,y3) != ccw(x1,y1, x2,y2,x4,y4);
        }
    }
}

