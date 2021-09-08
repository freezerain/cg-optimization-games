package simulacrum.marsLanderPuzzle;

class MarsLander {
    private final static double MARS_GRAVITY = 3.711; // m/ s^2
    GameState gs;
    public MarsLander(GameState gs) {
        this.gs = gs;
    }

    public void evaluateTurn(int rotateAngle, int power) {
        gs.power = Math.abs(power - gs.power) > 1 ? gs.power + (power - gs.power > 0 ? +1 : -1) :
                power;
        gs.rotateAngle = Math.abs(rotateAngle - gs.rotateAngle) > 15 ? gs.rotateAngle + (rotateAngle - gs.rotateAngle > 0 ? +15 : -15) :
                rotateAngle;
        double radians = Math.toRadians(gs.rotateAngle);;
        double xAcc = Math.sin(radians) * gs.power;
        double yAcc = Math.cos(radians) * gs.power - MARS_GRAVITY;
        gs.hSpeed -= xAcc;
        gs.vSpeed -= yAcc;
        gs.x = gs.x + gs.hSpeed - xAcc * 0.5;
        gs.y = gs.y - gs.vSpeed - yAcc * 0.5 - MARS_GRAVITY;
    }
}
