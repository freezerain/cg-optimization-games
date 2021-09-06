import java.util.*;

class Player {
    static final int[] FIBONACCI = {1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987,
            1597, 2584, 4181, 6765, 10946, 17711, 28657, 46368, 75025, 121393, 196418, 317811,
            514229, 832040, 1346269, 2178309, 3524578, 5702887, 9227465, 14930352, 24157817,
            39088169, 63245986, 102334155};

    public static void main(String[] args) {
        Scanner in = new Scanner(System.in);
        List<Chromosome> population = new ArrayList<>();
        Game game;
        GeneticAlgorithm ga = new GeneticAlgorithm();
        while (true) {
            Timer timer = new Timer();
            Actor player = new Actor(-1, in.nextInt(), in.nextInt());
            List<Actor> humans = new ArrayList<>();
            List<Actor> zombies = new ArrayList<>();
            int humanCount = in.nextInt();
            for (int i = 0; i < humanCount; i++)
                 humans.add(new Actor(in.nextInt(), in.nextInt(), in.nextInt()));
            int zombieCount = in.nextInt();
            for (int i = 0; i < zombieCount; i++){
                zombies.add(new Actor(in.nextInt(), in.nextInt(), in.nextInt()));
                in.nextInt();
                in.nextInt();
            }
           // printInputs(new Game(player, humans, zombies));
//            System.err.println("before eval:" + timer);
            ga.evaluate(new Game(player, humans, zombies), timer);
//            System.err.println("after eval:" + timer);
        }
    }

    static void printInputs(Game game) {
        System.err.println("---begin---");
        System.err.println(game.player.x + " " + game.player.y + " " + game.humans.size());
        game.humans.forEach(h -> System.err.println(h.id + " " + h.x + " " + h.y));
        System.err.println(game.zombies.size());
        game.zombies.forEach(z -> System.err.println(z.id + " " + z.x + " " + z.y + " 0 0 "));
        System.err.println("---end---");
    }

    static class GeneticAlgorithm {
        public static final int TOURNAMENT_SIZE = 5;
        private static final double ELITISM_PERCENTAGE = 20.0;
        private static final int PATH_LENGTH = 50;
        private static final int DESIRED_POPULATION_SIZE = 50;
        private static final double CROSSOVER_PERCENTAGE = 100.0 - ELITISM_PERCENTAGE;
        private static final double MUTATION_CHANCE = 2.0;
        private final long EVALUATE_TIME = 95;
        private final Random r = new Random();
        private List<Chromosome> pop = new ArrayList<>();

        public GeneticAlgorithm() {
            topUpWithEmptyChromosome(pop);
        }

        public void evaluate(Game game, Timer t) {
//            System.err.println("Population before top up:"+pop.size());
            if (pop.size() < DESIRED_POPULATION_SIZE) topUpWithEmptyChromosome(pop);
//            System.err.println("Population after top up:"+pop.size());
            int elitismIndex = (int) (pop.size() * (ELITISM_PERCENTAGE / 100));
            List<Chromosome> elite = pop.subList(0, elitismIndex);
//            System.err.println("elitism index:" + elitismIndex);
            List<List<Gene>> elitePath = new ArrayList<>();
            for (Chromosome c : elite){
                elitePath.add(c.path);
            }
            elite = simulateGenes(game, elitePath);
            //List<Chromosome> previousPopulation = pop.subList(elitismIndex, pop.size());
            //System.err.println("Elite size:" + elite.size() + " previous size:" + 
            // previousPopulation.size());
            int counter = 0;
            while (t.getTime() < EVALUATE_TIME) {
                //Crossover
                List<List<Gene>> crossoverGenes = crossoverPopulation(pop);
                //Mutate
                mutate(crossoverGenes);
                //Fitness
                pop = simulateGenes(game, crossoverGenes);
                pop.addAll(elite);
                counter++;
            }
            System.err.println("Genetic evaluations: " + counter);
//            System.err.println("Crossover population:" + pop.size());
            //Merge with elite
            //pop.addAll(elite);
//            System.err.println("population size:" + pop.size());
            //Sort
            pop.sort(Comparator.reverseOrder());
            //get best
            Chromosome bestChromosome = pop.get(0);
            //send message
            Gene bestGene = bestChromosome.path.get(0);
            int nextX = Math.max(0,
                    (int) (game.player.x + bestGene.radius * Math.cos(bestGene.angle * 22.5)));
            int nextY = Math.max(0,
                    (int) (game.player.y + bestGene.radius * Math.sin(bestGene.angle * 22.5)));
            System.err.println(bestChromosome);
            removeStep();
            if (bestChromosome.score == -1){
                Actor closestHuman = game.humans.get(0);
                int distance = Integer.MAX_VALUE;
                humanLoop:
                for (Actor h : game.humans){
                    int distanceToPlayer = (int) Math.sqrt(game.player.getDistanceSqrt(h));
                    System.err.println(distanceToPlayer);
                    distanceToPlayer-=2000;
                    int turnsToPlayer = distanceToPlayer/1000 + (distanceToPlayer%1000!=0? 1 : 0);
                    for (Actor z : game.zombies){
                        int distanceSqrt = (int) Math.sqrt(h.getDistanceSqrt(z));
                        int turnsToZombie = distanceSqrt / 400 + (distanceSqrt%400!=0? 1 : 0);
                        if(turnsToZombie<turnsToPlayer)continue humanLoop;
                    }
                    if(turnsToPlayer<distance){
                        distance = turnsToPlayer;
                        closestHuman = h;
                    }
                }
                System.err.println(distance);
                nextX = closestHuman.x;
                nextY = closestHuman.y;
            }
            System.out.println(nextX + " " + nextY);
            //remove step
            //save result
        }

        private void removeStep() {
            for (int i = pop.size() - 1; i >= 0; i--){
                Chromosome c = pop.get(i);
                c.path.remove(0);
                if (c.path.isEmpty()) pop.remove(i);
            }
        }

        private void topUpWithEmptyChromosome(List<Chromosome> population) {
            if (pop.size() >= DESIRED_POPULATION_SIZE) return;
            List<List<Gene>> randomGenes = getRandomGenes(DESIRED_POPULATION_SIZE - pop.size());
            randomGenes.forEach(l -> population.add(new Chromosome(l)));
        }

        private List<List<Gene>> getRandomGenes(int maxAmount) {
            List<List<Gene>> result = new ArrayList<>();
            for (int j = 0; j < maxAmount; j++){
                List<Gene> path = new ArrayList<>();
                for (int i = 0; i < PATH_LENGTH; i++){
                    double distanceWeight = r.nextDouble();
                    int angle = r.nextInt(16);
                    int distance = distanceWeight > 0.50 ? 1000 : (int) (2000 * distanceWeight);
                    Gene nextGene = new Gene(distance, angle);
                    path.add(nextGene);
                }
                result.add(path);
            }
            return result;
        }

        public List<Chromosome> simulateGenes(Game game, List<List<Gene>> path) {
            List<Chromosome> result = new ArrayList<>();
            for (List<Gene> genes : path){
                GameState gs = new GameState(game, genes);
                Chromosome chromosome = gs.evaluateGame();
                result.add(chromosome);
            }
            return result;
        }

        public void mutate(List<List<Gene>> population) {
            for (int i = 0; i < population.size(); i++){
                if (r.nextDouble() < MUTATION_CHANCE / 100){
                    int size = population.get(i).size();
                    List<Gene> mutatedGene = new ArrayList<>();
                    for (int j = 0; j < size; j++){
                        double distanceWeight = r.nextDouble();
                        int angle = r.nextInt(360);
                        int distance = distanceWeight < 0.25 ? 0 :
                                distanceWeight > 0.75 ? 1000 : (int) (1000 * distanceWeight);
                        Gene nextGene = new Gene(distance, angle);
                        mutatedGene.add(nextGene);
                    }
                    population.set(i, mutatedGene);
                }
            }
        }

        public List<List<Gene>> crossoverPopulation(List<Chromosome> population) {
            int maxAmount = (int) (population.size() * (CROSSOVER_PERCENTAGE / 100));
            List<List<Gene>> result = new ArrayList<>();
            for (int i = 0; i < maxAmount; i++){
                Chromosome p1 = selectTournament(population, TOURNAMENT_SIZE);
                Chromosome p2 = selectTournament(population, TOURNAMENT_SIZE);
                List<Gene> newGene = crossoverGene(p1.path, p2.path);
                result.add(newGene);
            }
            return result;
        }

        public List<Gene> crossoverGene(List<Gene> p1, List<Gene> p2) {
            if (p1.size() > p2.size()){
                List<Gene> temp = p1;
                p1 = p2;
                p2 = temp;
            }
            int crossoverIndex = r.nextInt(Math.min(p1.size(), p2.size()));
            List<Gene> child = new ArrayList<>();
            for (int i = 0; i < p2.size(); i++){
                child.add(i <= crossoverIndex ? p1.get(i) : p2.get(i));
            }
            return child;
        }

        public Chromosome selectTournament(List<Chromosome> population, int tournamentSize) {
            Chromosome bestInTournament = null;
            for (int j = 0; j < tournamentSize; j++){
                int index = r.nextInt(population.size());
                Chromosome chromosome = population.get(index);
                if (bestInTournament == null ||
                    bestInTournament.compareTo(chromosome) < 0) bestInTournament = chromosome;
            }
            return bestInTournament;
        }
    }

    static class Game {
        Actor player;
        List<Actor> humans;
        List<Actor> zombies;

        public Game(Actor player, List<Actor> humans, List<Actor> zombies) {
            this.player  = player;
            this.humans  = humans;
            this.zombies = zombies;
        }
    }

    static class GameState {
        Actor player;
        List<Actor> humans;
        List<Actor> zombies;
        int score = 0;
        List<Gene> path;
        ArrayList<ActorDistance> playerToZombies;
        ArrayList<ActorDistance> playerToHumans;
        HashMap<Actor, ArrayList<ActorDistance>> zombieToHumans;
        List<Actor> humanKillList;

        public GameState(Game game, List<Gene> path) {
            player = new Actor(game.player);
            humans = new ArrayList<>();
            game.humans.forEach(h -> humans.add(new Actor(h)));
            zombies = new ArrayList<>();
            game.zombies.forEach(h -> zombies.add(new Actor(h)));
            this.path = path;
        }

        public Chromosome evaluateGame() {
            for (Gene gene : path){
                //0) Precalculate distances
                //getAllDistances();
                //1) Move Zombies
                moveZombies();
                //2) Move Player
                movePlayer(gene);
                //3) Player kills zombies
                killZombies();
                //4) Zombies eat humans
                killHumans();
                if (zombies.isEmpty() || humans.isEmpty()) break;
            }
            return new Chromosome(path, score, humans.size(), zombies.size(), humansCanBeSaved());
        }

        /*private void getAllDistances() {
            playerToZombies = getSortedDistanceMaps(player, zombies, false);
            playerToHumans  = getSortedDistanceMaps(player, humans, false);
            zombieToHumans  = new HashMap<>();
            for (Actor z : zombies) zombieToHumans.put(z, getSortedDistanceMaps(z, humans, true));
        }

        private List<Actor> getZombieKillList() {
            ArrayList<Actor> result = new ArrayList<>();
            for (Actor z : zombies){
                if (z.getDistanceSqrt(player) <= 2000){
                    result.add(z);
                }
            }
            return result;
        }

        private ArrayList<ActorDistance> getSortedDistanceMaps(Actor caller, List<Actor> actors,
                                                               boolean addPlayer) {
            ActorDistance[] distances = new ActorDistance[actors.size() + (addPlayer ? 1 : 0)];

            ArrayList<ActorDistance> result = new ArrayList<>(actors.size() + (addPlayer ? 1 : 0));
            for (int i = 0; i < actors.size(); i++)
                 distances[i] = new ActorDistance(caller, actors.get(i),
                         actors.get(i).getDistanceSqrt(caller));
            if (addPlayer) distances[distances.length - 1] = new ActorDistance(caller, player,
                    player.getDistanceSqrt(caller));
            Arrays.sort(distances, Comparator.comparingInt(a -> a.distance));
            return new ArrayList<>(List.of(distances));
        }*/

        private int humansCanBeSaved() {
            return humans.size();
        }

        private void movePlayer(Gene move) {
            double vectorX = move.radius * Math.cos(move.angle * 22.5);
            double vectorY = move.radius * Math.sin(move.angle * 22.5);
            int nextX = Math.max(0, (int) (player.x + vectorX));
            int nextY = Math.max(0, (int) (player.y + vectorY));
            player.moveToNext(nextX, nextY, 1000);
        }

        private void killHumans() {
            for (Actor a : humanKillList){
                humans.remove(a);
            }
        }

        private void killZombies() {
            int fibSum = 0;
            int combo = 0;
            for (int i = zombies.size() - 1; i >= 0; i--){
                if (player.getDistanceSqrt(zombies.get(i)) <= 4000000){
                    fibSum += FIBONACCI[combo++];
                    zombies.remove(i);
                }
            }
            score += 10 * (humans.size() * humans.size()) * fibSum;
        }

        private void moveZombies() {
            humanKillList = new ArrayList<>();
            for (Actor z : zombies){
                Actor closest = player;
                int distance = player.getDistanceSqrt(z);
                for (Actor h : humans){
                    int hDis = z.getDistanceSqrt(h);
                    if (hDis < distance){
                        distance = hDis;
                        closest  = h;
                    }
                }
                if (closest.id != player.id && distance <= 160000) humanKillList.add(closest);
                z.moveToNext(closest.x, closest.y, 400);
            }
        }
    }

    static class ActorDistance {
        Actor from;
        Actor to;
        int distance;

        public ActorDistance(Actor from, Actor to, int distance) {
            this.from     = from;
            this.to       = to;
            this.distance = distance;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            ActorDistance that = (ActorDistance) o;
            if (from.equals(that.from) && to.equals(that.to)) return true;
            if (from.equals(that.to) && to.equals(that.from)) return true;
            return false;
        }

        @Override
        public int hashCode() {
            return 31 * (from.hashCode() + to.hashCode());
        }
    }

    static class Chromosome implements Comparable<Chromosome> {
        List<Gene> path;
        int score;
        int humansAlive;
        int zombieAlive;
        int fitnessScore;
        int humansCanBeSaved;

        public Chromosome(List<Gene> path) {
            this(path, 0, 0, 0, 0);
        }

        public Chromosome(List<Gene> path, int score, int humansAlive, int zombieAlive,
                          int humansCanBeSaved) {
            this.path             = path;
            this.score            = humansAlive == 0 ? -1 : score;
            this.humansAlive      = humansAlive;
            this.zombieAlive      = zombieAlive;
            this.humansCanBeSaved = humansCanBeSaved;
            calculateFitnessScore();
        }

        private void calculateFitnessScore() {
            if (humansAlive == 0){
                fitnessScore = -1;
                return;
            }
        }

        @Override
        public int compareTo(Chromosome o) {
            if (score != o.score) return Integer.compare(score, o.score);
            if (humansAlive != o.humansAlive) return Integer.compare(humansAlive, o.humansAlive);
            if (humansCanBeSaved != o.humansCanBeSaved)
                return Integer.compare(humansCanBeSaved, o.humansCanBeSaved);
            if (zombieAlive != o.zombieAlive) return Integer.compare(zombieAlive, o.zombieAlive);
            if (path.size() != o.path.size()) return Integer.compare(o.path.size(), path.size());
            return 0;
        }

        @Override
        public String toString() {
            return "Chromosome{" + "score=" + score + ", humansAlive=" + humansAlive +
                   ", zombieAlive=" + zombieAlive + ", humansCanBeSaved=" + humansCanBeSaved + '}';
        }
    }

    static class Gene {
        //0-1000
        int radius;
        //0-1000
        int angle;

        public Gene(int radius, int angle) {
            this.radius = radius;
            this.angle  = angle;
        }
    }

    static class Actor {
        int id;
        int x;
        int y;

        public Actor(int id, int x, int y) {
            this.id = id;
            this.x  = x;
            this.y  = y;
        }

        public Actor(Actor a) {
            this.id = a.id;
            this.x  = a.x;
            this.y  = a.y;
        }

        public int getDistanceSqrt(Actor target) {
            return (target.x - x) * (target.x - x) + (target.y - y) * (target.y - y);
        }

        public int getDistanceSqrt(int targetX, int targetY) {
            return (targetX - x) * (targetX - x) + (targetY - y) * (targetY - y);
        }

        public List<Actor> getInRange(List<Actor> actorList, int range) {
            List<Actor> result = new ArrayList<>();
            range = range * range;
            for (Actor a : actorList){
                if (this.getDistanceSqrt(a) <= range) result.add(a);
            }
            return result;
        }

        public Actor getClosest(List<Actor> actorList) {
            Actor closest = null;
            int distance = Integer.MAX_VALUE;
            for (Actor a : actorList){
                if (this.getDistanceSqrt(a) < distance) closest = a;
            }
            return closest;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            Actor actor = (Actor) o;
            return id == actor.id;
        }

        @Override
        public int hashCode() {
            return id;
        }

        public void moveToNext(int nextX, int nextY, int movementDistance) {
            double vectorX = nextX - x;
            double vectorY = nextY - y;
            double normaliser = Math.max(Math.abs(vectorX), Math.abs(vectorY));
            int newX = (int) (x + (vectorX / normaliser * movementDistance));
            int newY = (int) (y + (vectorY / normaliser * movementDistance));
            x = x < nextX ? Math.min(nextX, newX) : Math.max(nextX, newX);
            y = y < nextY ? Math.min(nextY, newY) : Math.max(nextY, newY);
        }
    }

    static class Timer {
        long start;

        public Timer() {
            start = System.currentTimeMillis();
        }

        public long getTime() {
            return System.currentTimeMillis() - start;
        }

        @Override
        public String toString() {
            return (System.currentTimeMillis() - start) + "ms.";
        }
    }
}