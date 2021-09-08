package simulacrum.marsLanderPuzzle;

class GameState {
    double x;
    double y;
    double hSpeed;
    double vSpeed;
    int fuel;
    int rotateAngle;
    int power;

    public GameState(double x, double y, double hSpeed, double vSpeed, int fuel, int rotateAngle,
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
}
