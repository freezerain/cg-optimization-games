package simulacrum.marsLanderPuzzle;

class MarsLander {
    private final static double MARS_GRAVITY = 3.711; // m/ s^2
    double x;
    double y;
    double hSpeed;
    double vSpeed;
    int fuel;
    int rotateAngle;
    int power;

    public MarsLander(double x, double y, double hSpeed, double vSpeed, int fuel, int rotateAngle,
                     int power) {
        this.x           = x;
        this.y           = y;
        this.hSpeed      = hSpeed;
        this.vSpeed      = vSpeed;
        this.fuel        = fuel;
        this.rotateAngle = rotateAngle;
        this.power       = power;
    }

    @Override
    public String toString() {
        return (int) Math.round(x) + " " + (int) Math.round(y) + " " + (int) Math.round(hSpeed) +
               " " + (int) Math.round(vSpeed) + " " + fuel + " " + rotateAngle + " " + power;
    }

    public void evaluateTurn(int desiredAngle, int desiredPower) {
        power = Math.abs(power - desiredPower) > 1 ? power + (desiredPower - power > 0 ? +1 : -1) :
                desiredPower;
        rotateAngle = Math.abs(rotateAngle - desiredAngle) > 15 ? rotateAngle + (desiredAngle - rotateAngle > 0 ? +15 : -15) :
                desiredAngle;
        double radians = Math.toRadians(rotateAngle);;
        double xAcc = Math.sin(radians) * power;
        double yAcc = Math.cos(radians) * power - MARS_GRAVITY;
        x += hSpeed - xAcc * 0.5;
        y += vSpeed + yAcc * 0.5;
        hSpeed -= xAcc;
        vSpeed += yAcc;
        
        fuel -= power;
    }
}
