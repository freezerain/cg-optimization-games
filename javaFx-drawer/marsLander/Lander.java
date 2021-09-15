package marsLander;

public class Lander {
    double x,y,hSpeed,vSpeed;
    int fuel, angle, power;
    boolean isLanded, isSafeLanded;

    public Lander(double x, double y, double hSpeed, double vSpeed, int fuel, int angle,
                  int power, boolean isLanded, boolean isSafeLanded) {
        this.x            = x;
        this.y            = y;
        this.hSpeed       = hSpeed;
        this.vSpeed       = vSpeed;
        this.fuel         = fuel;
        this.angle        = angle;
        this.power        = power;
        this.isLanded     = isLanded;
        this.isSafeLanded = isSafeLanded;
    }
}
