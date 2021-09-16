import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.util.*;

class Player {
    private final static int MARS_WIDTH = 7000;
    private final static int MARS_HEIGHT = 3000;
    public static GeneticAlgorithm GE = new GeneticAlgorithm();
    private static double[] landscape;
    private static double safeAreaBegin = -1;
    private static double safeAreaEnd = -1;
    private static int STARTING_FUEL = 3000;

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
            if (STARTING_FUEL == 3000) STARTING_FUEL = gameState.fuel;
            /*System.err.println(((int) gameState.x) + " " + ((int) gameState.y) + " " +
                               ((int) gameState.hSpeed) + " " + ((int) gameState.vSpeed) + " " +
                               gameState.fuel + " " + gameState.angle + " " + gameState.power);*/
            List<Individual> evaluate = GE.evaluate(gameState, t, new int[1]);
            if (!evaluate.isEmpty()){
                Individual individual = evaluate.get(0);
                GameState nextGS = individual.gameStateList.get(0);
                System.out.println((nextGS.angle) + " " + (nextGS.power));
            } else {
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

    static class GeneticAlgorithm {
        public static boolean IS_TOURNAMENT_SELECT = false;
        public static boolean IS_POINT_CROSSOVER = false;
        public static boolean RANDOM_CROSSOVER_ON_DUPLICATE = true;
        public static boolean REMOVE_DUPLICATES = true;
        public static boolean REMOVE_STEP = false;
        public static int TOURNAMENT_SIZE = 5;
        public static double ELITISM_PERCENTAGE = 0.2;
        public static int INDIVIDUAL_LENGTH = 200;
        public static int DESIRED_POPULATION_SIZE = 200;
        public static double CROSSOVER_PERCENTAGE = 0.6;
        public static double MUTATION_CHANCE = 0.02;
        public static double[] GENE_WEIGHTS = new double[]{0.3, 0.25, 0.25, 0.2, 1};
        private final Random r = new Random();
        public long EVALUATE_TIME = 100;
        private List<Individual> pop = new ArrayList<>();


        public List<Individual> evaluate(GameState game, long t, int[] counterArr) {
            if (pop.size() < DESIRED_POPULATION_SIZE) pop.addAll(
                    simulateGenes(game, getRandomGenes(DESIRED_POPULATION_SIZE - pop.size())));
            int counter = 0;
            while (System.currentTimeMillis() - t < EVALUATE_TIME) {
                List<List<Gene>> elite = new ArrayList<>();
                int elitismIndex = (int) (pop.size() * ELITISM_PERCENTAGE);
                for (int i = 0; i < elitismIndex; i++)
                     elite.add(pop.get(i).path);
                List<List<Gene>> newPop = crossover(pop);
                mutate(newPop);
                elite.addAll(newPop);
                newPop = elite;
                if (newPop.size() > DESIRED_POPULATION_SIZE) newPop = newPop.subList(0,
                        DESIRED_POPULATION_SIZE);
                else if (newPop.size() < DESIRED_POPULATION_SIZE) newPop.addAll(
                        getRandomGenes(DESIRED_POPULATION_SIZE - newPop.size()));
                pop = simulateGenes(game, newPop);
                if (REMOVE_DUPLICATES) pop = new ArrayList<>(new HashSet<>(pop));
                pop.sort(Comparator.comparingDouble((Individual i) -> i.fitnessScore).reversed());
                counter++;
            }
            counterArr[0] += counter;

            //System.err.println("Genetic evaluations: " + counter);
            //System.err.println(pop.get(0).fitnessScore);
           // System.err.println(pop.get(0).gameStateList.get(pop.get(0).gameStateList.size() - 1));
            if(REMOVE_STEP) removeStep();
            return pop;
        }

        //Rework remove step
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
            for (List<Gene> genes : population)
                for (int j = 0; j < genes.size(); j++)
                    if (r.nextDouble() < MUTATION_CHANCE) genes.set(j, new Gene());
        }

        public List<List<Gene>> crossover(List<Individual> population) {
            List<List<Gene>> result = new ArrayList<>();
            while (result.size() < DESIRED_POPULATION_SIZE * CROSSOVER_PERCENTAGE) {
                Individual ind1 = IS_TOURNAMENT_SELECT ? selectTournament(population) :
                        selectWheel(population);
                Individual ind2 = IS_TOURNAMENT_SELECT ? selectTournament(population) :
                        selectWheel(population);
                List<Gene> p1 = ind1.path;
                List<Gene> p2 = ind2.path;
                if (RANDOM_CROSSOVER_ON_DUPLICATE && p1 == p2){
                    p2 = new ArrayList<>();
                    for (int i = 0; i < p1.size(); i++)
                         p2.add(new Gene());
                } else if (p1 == p2){
                    do {
                        ind2 = selectTournament(population);
                    } while (p1 == p2);
                    p2 = ind2.path;
                }
                if (IS_POINT_CROSSOVER) crossoverOnPoint(result, p1, p2);
                else crossoverContinuous(result, p1, p2);
            }
            return result;
        }

        public void crossoverContinuous(List<List<Gene>> childList, List<Gene> p1, List<Gene> p2) {
            List<Gene> child1 = new ArrayList<>();
            List<Gene> child2 = new ArrayList<>();
            double crossoverPoint = r.nextDouble();
            double crossoverRest = 1.0 - crossoverPoint;
            for (int j = 0; j < p1.size() && j < p2.size(); j++){
                Gene parent1 = p1.get(j);
                Gene parent2 = p2.get(j);
                double angle1 = crossoverPoint * (parent1.angle + 15) +
                                crossoverRest * (parent2.angle + 15);
                double angle2 = crossoverRest * (parent1.angle + 15) +
                                crossoverPoint * (parent2.angle + 15);
                double power1 =
                        crossoverPoint * (parent1.power + 1) + crossoverRest * (parent2.power + 1);
                double power2 =
                        crossoverRest * (parent1.power + 1) + crossoverPoint * (parent2.power + 1);
                child1.add(new Gene((int) Math.round(angle1 - 15), (int) Math.round(power1 - 1)));
                child2.add(new Gene((int) Math.round(angle2 - 15), (int) Math.round(power2 - 1)));
            }
            childList.add(child1);
            childList.add(child2);
        }

        public void crossoverOnPoint(List<List<Gene>> childList, List<Gene> p1, List<Gene> p2) {
            int crossoverIndex = r.nextInt(Math.min(p1.size() - 2, p2.size() - 2)) + 1;
            List<Gene> child1 = new ArrayList<>(p1.subList(0, crossoverIndex));
            child1.addAll(p2.subList(crossoverIndex, p2.size()));
            List<Gene> child2 = new ArrayList<>(p2.subList(0, crossoverIndex));
            child2.addAll(p1.subList(crossoverIndex, p1.size()));
            childList.add(child1);
            childList.add(child2);
        }

        private Individual selectWheel(List<Individual> population) {
            double fitnessSum = 0.0;
            for (Individual i : population){
                fitnessSum += i.fitnessScore;
            }
            double[] probabilities = new double[population.size()];
            double previousProbability = 0.0;
            for (int i = 0; i < population.size(); i++){
                Individual ind = population.get(i);
                                   previousProbability += ind.fitnessScore / fitnessSum;
                probabilities[i] = previousProbability;
            }
            double selectProbability = r.nextDouble();
            for (int i = 0; i < probabilities.length; i++){
                if (selectProbability < probabilities[i]) return population.get(i);
            }
            return population.get(population.size() - 1);
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
            return fitnessScore - that.fitnessScore < 0.000001 &&
                   gameStateList.equals(that.gameStateList);
        }

        @Override
        public int hashCode() {
            int result;
            long temp;
            temp   = Double.doubleToLongBits(fitnessScore);
            result = (int) (temp ^ (temp >>> 32));
            result = 31 * result * gameStateList.hashCode();
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

            //Normalize to 0-100 range by dividing current / max
            double normDistance = 100.0 - distance / 7000.0 * 100.0;
            //Make score exponential
            normDistance *= normDistance;
            double normHS = 100.0 - (hSpeed < 20.0 ? 0 : hSpeed) / 500.0 * 100.0;
            normHS *= normHS;
            double normVS = 100.0 - (vSpeed < 40.0 ? 0 : vSpeed) / 500.0 * 100.0;
            normHS *= normVS;
            double normAngle = 100.0 - angle / 90.0 * 100.0;
            normAngle *= normAngle;
            double normFuel = fuel / (double) STARTING_FUEL * 100.0;
            normFuel *= normFuel;

            //Get fitness
            fitnessScore = normDistance * GeneticAlgorithm.GENE_WEIGHTS[0] +
                           normHS * GeneticAlgorithm.GENE_WEIGHTS[1] +
                           normVS * GeneticAlgorithm.GENE_WEIGHTS[2] +
                           normAngle * GeneticAlgorithm.GENE_WEIGHTS[3];
            //Add fuel to score only after safe landing
            if (lastState.isSafeLanded) fitnessScore += normFuel * GeneticAlgorithm.GENE_WEIGHTS[4];
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
            angle = angleDir == 0 ? -15 : angleDir == 1 ? 0 : 15;
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
            if (fuel <= 0 || isOutOfBounds() || isIntersect(oldX, oldY)){
                isLanded     = true;
                isSafeLanded = (!isOutOfBounds()) && fuel >= 0 && isSafeLand();
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

            if (Double.compare(gameState.x, x) != 0) return false;
            if (Double.compare(gameState.y, y) != 0) return false;
            if (Double.compare(gameState.hSpeed, hSpeed) != 0) return false;
            if (Double.compare(gameState.vSpeed, vSpeed) != 0) return false;
            if (fuel != gameState.fuel) return false;
            if (angle != gameState.angle) return false;
            if (power != gameState.power) return false;
            if (isLanded != gameState.isLanded) return false;
            return isSafeLanded == gameState.isSafeLanded;
        }

        @Override
        public int hashCode() {
            int result;
            long temp;
            temp   = Double.doubleToLongBits(x);
            result = (int) (temp ^ (temp >>> 32));
            temp   = Double.doubleToLongBits(y);
            result = 31 * result + (int) (temp ^ (temp >>> 32));
            temp   = Double.doubleToLongBits(hSpeed);
            result = 31 * result + (int) (temp ^ (temp >>> 32));
            temp   = Double.doubleToLongBits(vSpeed);
            result = 31 * result + (int) (temp ^ (temp >>> 32));
            result = 31 * result + fuel;
            result = 31 * result + angle;
            result = 31 * result + power;
            result = 31 * result + (isLanded ? 1 : 0);
            result = 31 * result + (isSafeLanded ? 1 : 0);
            return result;
        }

      /*  @Override
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
        }   */
    }
}

