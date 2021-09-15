import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.util.*;

class Player {
    private final static int MARS_WIDTH = 7000;
    private final static int MARS_HEIGHT = 3000;
    private static final GeneticAlgorithm GE = new GeneticAlgorithm();
    private static double[] landscape;
    private static double safeAreaBegin = -1;
    private static double safeAreaEnd = -1;

    public static void main(String args[]) {
        Scanner in = new Scanner(System.in);
        int N = in.nextInt(); // the number of points used to draw the surface of Mars.
        landscape = new double[N * 2];
        double lastX = -1;
        double lastY = -1;
        for (int i = 0; i < N * 2; i += 2){
            int x = in.nextInt();
            int y = in.nextInt();
            if (safeAreaBegin < 0 && y == lastY && x - lastX >= 1000){
                safeAreaBegin = lastX + 100;
                safeAreaEnd   = x - 100;
            }
            lastX            = x;
            lastY            = y;
            landscape[i]     = x;
            landscape[i + 1] = y;
        }

        for (double i : landscape){
            System.err.print(i + " ");
        }
        System.err.println();

        while (true) {
            long t = System.currentTimeMillis();
            GameState gameState = new GameState(in.nextInt(), in.nextInt(), in.nextInt(),
                    in.nextInt(), in.nextInt(), in.nextInt(), in.nextInt());
            /*System.err.println(((int) gameState.x) + " " + ((int) gameState.y) + " " +
                               ((int) gameState.hSpeed) + " " + ((int) gameState.vSpeed) + " " +
                               gameState.fuel + " " + gameState.angle + " " + gameState.power);*/
            List<Individual> evaluate = GE.evaluate(gameState, t);
            if(!evaluate.isEmpty()) {
                Individual individual = evaluate.get(0);
                GameState nextGS = individual.gameStateList.get(0);
                System.out.println((nextGS.angle) + " " + (nextGS.power));
            }else {
                System.out.println(0 + " " + 0);
            }
            //Gene gene = individual.path.get(0);
           // System.out.println(-50 + " " + 1);
        }
    }

    static double[] getNearestLandscapeSegment(double x) {
        for (int i = 2; i < landscape.length; i += 2){
            if (landscape[i] >= x)
                return new double[]{landscape[i - 2], landscape[i - 1], landscape[i],
                        landscape[i + 1]};
        }
        return new double[]{0.0, 0.0, 0.0, 0.0};
    }

    void simulationInit(double[] landscape) {
        Player.landscape = landscape;
        double lastX = -1;
        double lastY = -1;
        for (int i = 0; i < landscape.length; i += 2){
            if (safeAreaBegin < 0 && landscape[i + 1] == lastY && landscape[i] - lastX >= 1000){
                safeAreaBegin = lastX + 100;
                safeAreaEnd   = landscape[i] - 100;
            }
            lastX = landscape[i];
            lastY = landscape[i + 1];
        }
    }

    List<Individual> simulationStep(GameState gs) {
        return GE.evaluate(gs, System.currentTimeMillis());
    }

    static class GeneticAlgorithm {
        public static final int TOURNAMENT_SIZE = 5;
        private static final double ELITISM_PERCENTAGE = 10.0;
        private static final int INDIVIDUAL_LENGTH = 200;
        private static final int DESIRED_POPULATION_SIZE = 100;
        private static final double CROSSOVER_PERCENTAGE = 70.0;
        private static final double MUTATION_CHANCE = 2.0;
        private final long EVALUATE_TIME = 95;
        private final Random r = new Random();
        private List<Individual> pop = new ArrayList<>();

        public List<Individual> evaluate(GameState game, long t) {
            if (pop.size() < DESIRED_POPULATION_SIZE) pop.addAll(
                    simulateGenes(game, getRandomGenes(DESIRED_POPULATION_SIZE - pop.size())));
            int elitismIndex = (int) (pop.size() * (ELITISM_PERCENTAGE / 100));
            List<List<Gene>> elitePath = new ArrayList<>();
            for (int i = 0; i < elitismIndex; i++)
                 elitePath.add(pop.get(i).path);
            List<Individual> elite = simulateGenes(game, elitePath);
            
            int counter = 0;
            while (System.currentTimeMillis() - t < EVALUATE_TIME) {
                elitePath = new ArrayList<>();
                elitismIndex = (int) (pop.size() * (ELITISM_PERCENTAGE / 100));
                for (int i = 0; i < elitismIndex; i++)
                     elitePath.add(pop.get(i).path);
                List<List<Gene>> crossoverGenes = crossoverPopulation(pop);
                int newPopSize = crossoverGenes.size() + elite.size();
                
                crossoverGenes.addAll(elitePath);
                mutate(crossoverGenes);
                if (newPopSize > DESIRED_POPULATION_SIZE) crossoverGenes = crossoverGenes.subList(0,
                        DESIRED_POPULATION_SIZE);
                else
                if (newPopSize < DESIRED_POPULATION_SIZE) crossoverGenes.addAll(
                        getRandomGenes(DESIRED_POPULATION_SIZE - newPopSize));
                pop = simulateGenes(game, crossoverGenes);
                pop = new ArrayList<>(new HashSet<>(pop));
                pop.sort(Comparator.comparingDouble((Individual i) -> i.fitnessScore).reversed());
                //if(pop.size()>DESIRED_POPULATION_SIZE - elitismIndex) pop = pop.subList(0, DESIRED_POPULATION_SIZE-elitismIndex); 
                counter++;
            }


           // pop.addAll(elite);
            System.err.println("Genetic evaluations: " + counter);
          //  pop.sort(Comparator.comparingDouble((Individual i) -> i.fitnessScore).reversed());
            System.err.println(pop.get(0).fitnessScore);
            System.err.println(pop.get(0).gameStateList.get(pop.get(0).gameStateList.size() - 1));
            pop.get(0).getFitness();
            removeStep();
            return pop;
        }

        private void removeStep() {
            for (int i = pop.size() - 1; i >= 0; i--){
                Individual c = pop.get(i);
                if (!c.path.isEmpty()) c.path.remove(0);
                if (!c.gameStateList.isEmpty()) c.gameStateList.remove(0);
                if (c.path.isEmpty() || c.gameStateList.isEmpty()) pop.remove(i);
            }
        }

        private List<List<Gene>> getRandomGenes(int maxAmount) {
            List<List<Gene>> result = new ArrayList<>();
            for (int j = 0; j < maxAmount; j++){
                List<Gene> path = new ArrayList<>();
                for (int i = 0; i < INDIVIDUAL_LENGTH; i++)
                     path.add(new Gene());
                result.add(path);
            }
            return result;
        }

        public List<Individual> simulateGenes(GameState gs, List<List<Gene>> path) {
            List<Individual> result = new ArrayList<>();
            for (List<Gene> genes : path){
                ArrayList<GameState> history = new ArrayList<>();
                GameState newGS = new GameState(gs);
                for (Gene gene : genes){
                    newGS.simulate(gene.angle, gene.power);
                    history.add(new GameState(newGS));
                    if (newGS.isLanded) break;
                }
                while (history.size() < INDIVIDUAL_LENGTH && !newGS.isLanded) {
                    Gene gene = new Gene();
                    newGS.simulate(gene.angle, gene.power);
                    genes.add(gene);
                    history.add(new GameState(newGS));
                }
                result.add(new Individual(genes.subList(0, history.size()), history));
            }
            return result;
        }

        public void mutate(List<List<Gene>> population) {
            for (int i = 0; i < population.size(); i++){
                if (r.nextDouble() < MUTATION_CHANCE / 100){
                    int size = population.get(i).size();
                    int mutatedGeneIndex = r.nextInt(size);
                    population.get(i).set(mutatedGeneIndex, new Gene());
                }
            }
        }

        public List<List<Gene>> crossoverPopulation(List<Individual> population) {
            List<List<Gene>> result = new ArrayList<>();
            while (result.size() < DESIRED_POPULATION_SIZE * (CROSSOVER_PERCENTAGE / 100.0)) {
                Individual p1 = selectTournament(population);
                Individual p2 = null;
                do {
                    p2 = selectTournament(population);
                } while (p1 == p2);

                List<Gene> newGene = crossoverGene(p1.path, p2.path);
                result.add(newGene);

                //Continuous Genetic Algorithm
           /*     List<Gene> child1 = new ArrayList<>();
                List<Gene> child2 = new ArrayList<>();
                double crossoverPoint = r.nextDouble();
                double crossoverRest = 1.0 - crossoverPoint;
                for (int j = 0; j < p1.path.size() && j < p2.path.size(); j++){
                    Gene parent1 = p1.path.get(j);
                    Gene parent2 = p2.path.get(j);
                    double angle1 = crossoverPoint * (parent1.angle + 15) +
                                    crossoverRest * (parent2.angle + 15);
                    double angle2 = crossoverRest * (parent1.angle + 15) +
                                    crossoverPoint * (parent2.angle + 15);
                    double power1 = crossoverPoint * parent1.power + crossoverRest * parent2.power;
                    double power2 = crossoverRest * parent1.power + crossoverPoint * parent2.power;
                    child1.add(new Gene((int) Math.round(angle1 - 15), (int) Math.round(power1)));
                    child2.add(new Gene((int) Math.round(angle2 - 15), (int) Math.round(power2)));
                }
                result.add(child1);
                result.add(child2);*/

            }
            return result;
        }

        public List<Gene> crossoverGene(List<Gene> p1, List<Gene> p2) {
            if (p1.size() > p2.size()){
                List<Gene> temp = p1;
                p1 = p2;
                p2 = temp;
            }
            int crossoverIndex = r.nextInt(p1.size());
            List<Gene> child = new ArrayList<>();
            for (int i = 0; i < p2.size(); i++){
                child.add(i <= crossoverIndex ? p1.get(i) : p2.get(i));
            }
            return child;
        }

        public Individual selectTournament(List<Individual> population) {
            Individual bestInTournament = null;
            for (int j = 0; j < TOURNAMENT_SIZE; j++){
                int index = r.nextInt(population.size());
                Individual ind = population.get(index);
                if (bestInTournament == null ||
                    ind.fitnessScore > bestInTournament.fitnessScore) bestInTournament = ind;
            }
            return bestInTournament;
        }
    }

    static class Individual {
        List<Gene> path;
        List<GameState> gameStateList;
        double fitnessScore;

        public Individual(List<Gene> path, List<GameState> gameStateList) {
            this.path          = path;
            this.gameStateList = gameStateList;
            getFitness();
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            Individual that = (Individual) o;
            return gameStateList.equals(that.gameStateList);
        }

        @Override
        public int hashCode() {
            int result;
            result = 31 * gameStateList.hashCode();
            return result;
        }

        public void getFitness() {
            GameState lastState = gameStateList.get(gameStateList.size() - 1);
            fitnessScore = 0;
            double distance = lastState.x < safeAreaBegin ? safeAreaBegin - lastState.x :
                    lastState.x > safeAreaEnd ? lastState.x - safeAreaEnd : 0.0;
            double hSpeed = Math.abs(lastState.hSpeed);
            double vSpeed = Math.abs(lastState.vSpeed);
            int angle = Math.abs(lastState.angle);
            int fuel = Math.abs(lastState.fuel);

/*            //Distance 0-7000
            fitnessScore += 1.0 - (distance==0.0? 0 : 0.5 + 0.5/7000.0*distance);
            //speed -500 - 500
            fitnessScore += 1.0 - (hSpeed<=20.0? 0 : 0.5 + 0.5/500.0*hSpeed);
            fitnessScore += 1.0 - (vSpeed<=40.0? 0 : 0.5 + 0.5/500.0*vSpeed);
            //angle -90 - 90
            fitnessScore += 1.0 - (angle==0? 0 : 0.5 + 0.5/90.0*angle);
            //fuel 0 - 2000*/

  /*          //Distance 0-7000
            fitnessScore = 10000;
            fitnessScore -= distance;
            //speed -500 - 500
            fitnessScore -= hSpeed<=20.0? 0 : hSpeed ;
            fitnessScore -= vSpeed<=40.0? 0 : vSpeed ;
            //angle -90 - 90
            fitnessScore -= angle==0? 0 : angle ;
            
           // fitnessScore = Math.pow(fitnessScore, 2);
            
            */

            //Distance 0-7000
            fitnessScore = 100000000;
            fitnessScore -= Math.pow(distance,2);
            //speed -500 - 500
            fitnessScore -= Math.pow(hSpeed<=15.0? 0 : hSpeed ,3);
            fitnessScore -= Math.pow(vSpeed<=35.0? 0 : vSpeed, 3);
            //angle -90 - 90
            fitnessScore -= Math.pow(angle,1);
            if(!lastState.isSafeLanded) fitnessScore -= 10000000;
            
            fitnessScore += Math.pow(fuel,0.5); 
            
            
            // fitnessScore = Math.pow(fitnessScore, 2);
            
            
            
            
            
            
            
/*            if(lastState.isSafeLanded){
                fitnessScore= 4 + 1-(lastState.fuel/2000.0);
            }else{
                fitnessScore += 1.0 - ((lastState.x < safeAreaBegin ? safeAreaBegin - lastState.x :
                        lastState.x > safeAreaEnd ? lastState.x - safeAreaEnd : 0) / 6000.0);
                //speed -500 - 500
                fitnessScore += 1.0 - ((Math.abs(lastState.hSpeed)*Math.abs(lastState.hSpeed) / (500.0*500.0)));
                fitnessScore += 1.0 - ((Math.abs(lastState.vSpeed)*Math.abs(lastState.vSpeed) / (500.0*500.0)));
                //angle -90 - 90
                fitnessScore += 1 - (lastState.angle == 0 ? 0 : (Math.abs(lastState.angle) / 90.0));
            }*/
            
            
            
            /*fitnessScore = 1000;
            fitnessScore -= lastState.x <safeAreaBegin? safeAreaBegin - lastState.x : lastState.x>safeAreaEnd? lastState.x - safeAreaEnd : 0;
            fitnessScore -= Math.abs(lastState.hSpeed);
            fitnessScore -= Math.abs(lastState.vSpeed);
            fitnessScore -= Math.abs(lastState.angle);*/
            
            
            
           /* //Distance 0-7000
            fitnessScore += 1.0 - ((lastState.x < safeAreaBegin ? safeAreaBegin - lastState.x :
                    lastState.x > safeAreaEnd ? lastState.x - safeAreaEnd : 0) / 7000);
            //speed -500 - 500
            fitnessScore += 1.0 - (Math.abs(lastState.hSpeed) <= 20 ? 0 : (Math.abs(lastState
            .hSpeed) / 500.0));
            fitnessScore += 1.0 - (Math.abs(lastState.vSpeed) <= 40 ? 0 : (Math.abs(lastState
            .vSpeed) / 500.0));
            //angle -90 - 90
            fitnessScore += 1 - (lastState.angle == 0 ? 0 : (Math.abs(lastState.angle) / 90.0));
            //fuel 0 - 2000
            fitnessScore += lastState.fuel / 2000.0;
*/

            //Distance 0-7000
         /*   fitnessScore += (7000.0 - (lastState.x < safeAreaBegin ? safeAreaBegin - lastState.x :
                    lastState.x > safeAreaEnd ? lastState.x - safeAreaEnd : 0)) / 7000.0;
            //speed -500 - 500
            fitnessScore *= ((Math.abs(lastState.hSpeed) <= 20.0 ? 1.0 :
                    ((500.0 - Math.abs(lastState.hSpeed)) / 500.0)) +
                             (Math.abs(lastState.vSpeed) <= 40.0 ? 1.0 :
                                     ((500.0 - Math.abs(lastState.vSpeed)) / 500.0))) / 2.0;
            //angle -90 - 90
            fitnessScore *= (90.0 - Math.abs(lastState.angle)) / 90.0;
            //fuel 0 - 2000
            fitnessScore *= (2000.0 - lastState.fuel) / 2000.0;*/

        }
    }

    static class Gene {
        final public int angle;
        final public int power;

        public Gene(int angle, int power) {
            this.angle = angle;
            this.power = power;
        }

        public Gene() {
            Random r = new Random();
            int angleDir = r.nextInt(3);
            angle = angleDir==0? -15 : angleDir ==1? 0 : 15;
            power = r.nextInt(3) - 1;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            Gene gene = (Gene) o;

            if (angle != gene.angle) return false;
            return power == gene.power;
        }

        @Override
        public int hashCode() {
            int result = angle;
            result = 31 * result + power;
            return result;
        }
    }

    static class GameState {
        double x;
        double y;
        double hSpeed;
        double vSpeed;
        int fuel;
        int angle;
        int power;
        boolean isLanded = false;
        boolean isSafeLanded = false;

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
            x            = gs.x;
            y            = gs.y;
            hSpeed       = gs.hSpeed;
            vSpeed       = gs.vSpeed;
            fuel         = gs.fuel;
            angle        = gs.angle;
            power        = gs.power;
            isLanded     = gs.isLanded;
            isSafeLanded = gs.isSafeLanded;
        }

        public void simulate(int dAngle, int dPower) {
            if (dPower != 0) power = Math.max(0, Math.min(4, power + dPower));
            if (dAngle != 0) angle = Math.max(-90, Math.min(90, angle + dAngle));
            double oldX = x;
            double oldY = y;
            move();
            if (fuel<= 0 ||isOutOfBounds() || isIntersect(oldX, oldY)){
                isLanded     = true;
                isSafeLanded = (!isOutOfBounds()) && fuel >= 0 && isSafeLand() ;
            }
        }

        private boolean isOutOfBounds() {
            return x < 0.0 || x > 7000.0 || y < 0.0 || y > 3000.0;
        }

        private boolean isIntersect(double oldX, double oldY) {
            Line2D line = new Line2D.Double(new Point2D.Double(oldX, oldY),
                    new Point2D.Double(x, y));
            for (int i = 0; i < landscape.length - 2; i += 2){
                if (line.intersectsLine(landscape[i], landscape[i + 1], landscape[i + 2],
                        landscape[i + 3])) return true;
            }
            return false;
        }

        private boolean isSafeLand() {
            if (x > safeAreaEnd || x < safeAreaBegin || angle != 0 || Math.abs(vSpeed) > 40 ||
                Math.abs(hSpeed) > 20) return false;
            return true;
        }

        private void move() {
            double radians = Math.toRadians(angle);
            double xAcc = Math.sin(radians) * power;
            double yAcc = (Math.cos(radians) * power - 3.711);
            x += hSpeed - xAcc * 0.5;
            y += vSpeed + yAcc * 0.5;
            hSpeed -= xAcc;
            vSpeed += yAcc;
            fuel -= power;
        }

        @Override
        public String toString() {
            return "GameState{" + "x=" + x + ", y=" + y + ", hSpeed=" + hSpeed + ", vSpeed=" +
                   vSpeed + ", fuel=" + fuel + ", angle=" + angle + ", power=" + power +
                   ", isLanded=" + isLanded + ", isSafeLanded=" + isSafeLanded + '}';
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            GameState gameState = (GameState) o;
            if (isLanded != gameState.isLanded) return false;
            if (isSafeLanded != gameState.isSafeLanded) return false;
            if ((int) gameState.x != (int) x) return false;
            if ((int) gameState.y != (int) y) return false;
            if ((int) gameState.hSpeed != (int) hSpeed) return false;
            if ((int) gameState.vSpeed != (int) vSpeed) return false;
            if (fuel != gameState.fuel) return false;
            if (angle != gameState.angle) return false;

            return power == gameState.power;
        }

        @Override
        public int hashCode() {
            int result;
            result = (int) x;
            result = 31 * result + (int) y;
            result = 31 * result + (int) hSpeed;
            result = 31 * result + (int) vSpeed;
            result = 31 * result + fuel;
            result = 31 * result + angle;
            result = 31 * result + power;
            result = 31 * result + (isLanded ? 1 : 0);
            result = 31 * result + (isSafeLanded ? 1 : 0);
            return result;
        }
    }
}

