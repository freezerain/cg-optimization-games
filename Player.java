import java.util.*;

class Player {
    static final int[] FIBONACCI = {1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987,
            1597, 2584, 4181, 6765, 10946, 17711, 28657, 46368, 75025, 121393, 196418, 317811,
            514229, 832040, 1346269, 2178309, 3524578, 5702887, 9227465, 14930352, 24157817,
            39088169, 63245986, 102334155};

    public static void main(String[] args) {
        Scanner in = new Scanner(System.in);
        Genetic.Settings settings = new Genetic.Settings();
        Genetic.Solver solver = new Genetic.Solver(settings);
        List<Genetic.Individual> pop = new ArrayList<>();
        while (true) {
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
/*            System.err.println("Player: " + player.x + ", " + player.y);
            System.err.print("Humans: ");
            for (Actor human : humans){
                System.err.print(human.x+","+human.y+",");
            }
            System.err.println();
            System.err.print("Zombies: ");
            for (Actor human : zombies){
                System.err.print(human.x+","+human.y+",");
            }
            System.err.println();*/

            int[] counter = new int[1];
            Genetic.State startState = new Genetic.State(player, humans, zombies);
            pop = solver.evolve(pop, startState, System.currentTimeMillis(), counter);
            System.err.println("crossovers: " + counter[0]);
            Genetic.Individual bestInd = pop.get(0);
            System.err.println("ind: " + bestInd.fitness);
            int nextX;
            int nextY;
            if (bestInd.state.humans.isEmpty()){
                System.err.println("Solution not found! Simple logic.");
                Actor closestHuman = startState.humans.get(0);
                int distance = Integer.MAX_VALUE;
                humanLoop:
                for (Actor h : startState.humans){
                    int distanceToPlayer = (int) Math.sqrt(startState.player.getDistanceSqrt(h));
                    distanceToPlayer -= 2000;
                    int turnsToPlayer =
                            distanceToPlayer / 1000 + (distanceToPlayer % 1000 != 0 ? 1 : 0);
                    for (Actor z : startState.zombies){
                        int distanceSqrt = (int) Math.sqrt(h.getDistanceSqrt(z));
                        int turnsToZombie = distanceSqrt / 400 + (distanceSqrt % 400 != 0 ? 1 : 0);
                        if (turnsToZombie < turnsToPlayer) continue humanLoop;
                    }
                    if (turnsToPlayer < distance){
                        distance     = turnsToPlayer;
                        closestHuman = h;
                    }
                }
                nextX = closestHuman.x;
                nextY = closestHuman.y;
            } else {
                Genetic.Gene nextGene = bestInd.getGenes().get(0);
                Genetic.State nextState = bestInd.states.get(0);
                nextX = Math.max(0, (int) (startState.player.x + nextGene.getValues()[0] * Math.cos(
                        nextGene.getValues()[1] * 22.5)));
                nextY = Math.max(0, (int) (startState.player.y + nextGene.getValues()[0] * Math.sin(
                        nextGene.getValues()[1] * 22.5)));
            }
            solver.removeFirstStep(pop);
            System.out.println(nextX + " " + nextY);
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

        public void moveToNext(int nextX, int nextY) {
            x = nextX;
            y = nextY;
        }

        public void moveToNext(int nextX, int nextY, int movementDistance) {
            double vectorX = nextX - x;
            double vectorY = nextY - y;
            //double vectorX = nextX - x;
            //double vectorY = nextY - y;
            double dist = Math.sqrt(vectorX * vectorX + vectorY * vectorY);
            vectorX /= dist;
            vectorY /= dist;
            x += vectorX * movementDistance;
            y += vectorY * movementDistance;
            //double normaliser = Math.max(Math.abs(vectorX), Math.abs(vectorY));
            //int newX = (int) (x + (vectorX / normaliser * movementDistance));
            //int newY = (int) (y + (vectorY / normaliser * movementDistance));
            //x = x < nextX ? Math.min(nextX, newX) : Math.max(nextX, newX);
            //y = y < nextY ? Math.min(nextY, newY) : Math.max(nextY, newY);

        }
    }

    static class Genetic {
        private static final Random R = new Random();

        enum CrossoverType {
            CONTINUOUS {
                @Override
                public void doCrossover(List<Individual> childList, Individual p1, Individual p2) {
                    double crossoverPoint = R.nextDouble();
                    double crossoverRest = 1.0 - crossoverPoint;
                    List<Gene> genes1 = p1.getGenes();
                    List<Gene> genes2 = p2.getGenes();
                    List<Gene> newGenes1 = new ArrayList<>(genes1.size());
                    List<Gene> newGenes2 = new ArrayList<>(genes2.size());
                    int lastCrossover = Math.min(genes1.size(), genes2.size());
                    for (int j = 0; j < lastCrossover; j++){
                        Gene gene1 = genes1.get(j);
                        Gene gene2 = genes2.get(j);
                        int[] vars1 = gene1.getValues();
                        int[] vars2 = gene2.getValues();
                        int[] newVars1 = new int[vars1.length];
                        int[] newVars2 = new int[vars1.length];
                        for (int i = 0; i < vars1.length; i++){
                            newVars1[i] = (int) (crossoverPoint * vars1[i] +
                                                 crossoverRest * vars2[i]);
                            newVars2[i] = (int) (crossoverRest * vars1[i] +
                                                 crossoverPoint * vars2[i]);
                        }
                        newGenes1.add(new Gene(newVars1));
                        newGenes2.add(new Gene(newVars2));
                    }
                    childList.add(new Individual(newGenes1));
                    childList.add(new Individual(newGenes2));
                }
            }, POINT {
                @Override
                public void doCrossover(List<Individual> childList, Individual p1, Individual p2) {
                    List<Gene> genes1 = p1.getGenes();
                    List<Gene> genes2 = p2.getGenes();
                    if (genes1.size() < 2 || genes2.size() < 2){
                        childList.add(p1);
                        childList.add(p2);
                        return;
                    }
                    int crossoverIndex = (int) (R.nextDouble() *
                                                Math.min(genes1.size(), genes2.size()));
                    int max = Math.max(genes1.size(), genes2.size());
                    List<Gene> newGenes1 = new ArrayList<>(genes1.subList(0, crossoverIndex));
                    newGenes1.addAll(
                            new ArrayList<>(genes2.subList(crossoverIndex, genes2.size())));
                    List<Gene> newGenes2 = new ArrayList<>(genes2.subList(0, crossoverIndex));
                    newGenes1.addAll(
                            new ArrayList<>(genes1.subList(crossoverIndex, genes1.size())));
                    childList.add(new Individual(newGenes1));
                    childList.add(new Individual(newGenes2));
                }
            };

            public abstract void doCrossover(List<Individual> childList, Individual p1,
                                             Individual p2);
        }

        enum SelectType {
            TOURNAMENT {
                @Override
                public Individual doSelect(List<Individual> pop, int var) {
                    Individual bestInTournament = null;
                    for (int j = 0; j < var; j++){
                        int index = R.nextInt(pop.size());
                        Individual ind = pop.get(index);
                        if (bestInTournament == null || ind.getFitness() >
                                                        bestInTournament.getFitness()) bestInTournament = ind;
                    }
                    return bestInTournament;
                }
            }, WHEEL {
                @Override
                public Individual doSelect(List<Individual> pop, int var) {
                    double fitnessSum = 0.0;
                    for (Individual i : pop) fitnessSum += i.getFitness();
                    double[] probabilities = new double[pop.size()];
                    double previousProbability = 0.0;
                    for (int i = 0; i < pop.size(); i++){
                        Individual ind = pop.get(i);
                                           previousProbability += ind.getFitness() / fitnessSum;
                        probabilities[i] = previousProbability;
                    }
                    double selectProbability = Genetic.R.nextDouble();
                    for (int i = 0; i < probabilities.length; i++)
                        if (selectProbability < probabilities[i]) return pop.get(i);
                    return pop.get(0);
                }
            };

            public abstract Individual doSelect(List<Individual> pop, int var);
        }

        public static class State {
            Actor player;
            List<Actor> humans;
            List<Actor> zombies;
            long score = 0;
            KdTree.Node root;

            public State(Actor player, List<Actor> humans, List<Actor> zombies) {
                this.player = new Actor(player);
                this.humans = new ArrayList<>(humans.size());
                for (Actor h : humans){
                    this.humans.add(new Actor(h));
                }
                this.zombies = new ArrayList<>(zombies.size());
                for (Actor z : zombies){
                    this.zombies.add(new Actor(z));
                }
            }

            public State(State gs) {
                this(gs.player, gs.humans, gs.zombies);
                this.score = gs.score;
            }

            public void simulate(Gene gene) {
                // List<Actor> humanKillList = moveZombies();
                buildKdTree();
                moveZombies();
                movePlayer(gene.getValues()[0], gene.getValues()[1]);
                killZombies();
                // killHumans(humanKillList);
                killHumans();
            }
            
            private void buildKdTree(){
                root = null;
                root = KdTree.insert(root, player, true);
                for (Actor h : humans){
                    root = KdTree.insert(root, h, true);
                }
             /*   for (Actor z : zombies){
                    root = KdTree.insert(root, z, false);
                }*/
            }

            private void movePlayer(int radius, int angle) {
                double vectorX = radius * Math.cos(angle * 22.5);
                double vectorY = radius * Math.sin(angle * 22.5);
                int nextX = Math.min(15999, Math.max(0, (int) (player.x + vectorX)));
                int nextY = Math.min(8999, Math.max(0, (int) (player.y + vectorY)));
                player.moveToNext(nextX, nextY);
            }

            private void killHumans(/*List<Actor> humanKillList*/) {
                for (int i = humans.size() - 1; i >= 0; i--){
                    Actor h = humans.get(i);
                    for (Actor z : zombies){
                        if (h.x == z.x && h.y == z.y){
                            humans.remove(i);
                            break;
                        }
                    }
                }
            /*    for (Actor a : humanKillList)
                    humans.remove(a);*/
            }

            private void killZombies() {
                int fibSum = 0;
                int combo = 0;
                for (int i = zombies.size() - 1; i >= 0; i--){
                    if (player.getDistanceSqrt(zombies.get(i)) <= 4000000){
                        if (combo < 39) fibSum += FIBONACCI[combo++];
                        zombies.remove(i);
                    }
                }
                score += 10 * (humans.size() * humans.size()) * fibSum;
            }

            private void moveZombies() {
                //ArrayList<Actor> humanKillList = new ArrayList<>();
                for (Actor z : zombies){
                    Actor closest = KdTree.nearestNeighbour(root, z, true).actor;
                    //Actor closest = player;
                    /*int distance = player.getDistanceSqrt(z);s
                    for (Actor h : humans){
                        int hDis = z.getDistanceSqrt(h);
                        if (hDis < distance){
                            distance = hDis;
                            closest  = h;
                        }
                    }*/
                    //if (closest.id != player.id && distance <= 160000) humanKillList.add(closest);
                    z.moveToNext(closest.x, closest.y, 400);
                }
                //return humanKillList;
            }

            @Override
            public boolean equals(Object o) {
                if (this == o) return true;
                if (o == null || getClass() != o.getClass()) return false;

                State state = (State) o;

                if (score != state.score) return false;
                if (!player.equals(state.player)) return false;
                if (!humans.equals(state.humans)) return false;
                return zombies.equals(state.zombies);
            }

            @Override
            public int hashCode() {
                int result = player.hashCode();
                result = 31 * result + humans.hashCode();
                result = 31 * result + zombies.hashCode();
                result = 31 * result + (int) (score ^ (score >>> 32));
                return result;
            }
        }

        public static class Individual {
            List<Gene> genes = new ArrayList<>();
            List<State> states = new ArrayList<>();
            State state;
            double fitness;
            int zombiesMax;
            int humansMax;

            public Individual(List<Gene> genes) {
                this.genes = genes;
            }

            public Individual() {

            }

            public void removeFirstStep() {
                if (!genes.isEmpty()){
                    genes.remove(0);
                }
            }

            List<Gene> getGenes() {
                return genes;
            }

            double getFitness() {
                return fitness;
            }

            public void simulate(State state, Settings settings) {
                state      = new State(state);
                zombiesMax = state.zombies.size();
                humansMax  = state.humans.size();
                int i = 0;
                while (i <= settings.INDIVIDUAL_LENGTH && (!state.humans.isEmpty()) &&
                       (!state.zombies.isEmpty())) {
                    if (i >= this.genes.size()) genes.add(new Gene());
                    ;
                    state.simulate(this.genes.get(i));
                    i++;
                }
                this.state = state;
                fitness    = calculateFitness(settings);
            }

            private double calculateFitness(Settings settings) {
                int fibSum = 0;
                int combo = 0;
                for (int i = zombiesMax - 1; i >= 0; i--){
                    if (combo < 39) fibSum += FIBONACCI[combo++];
                }
                long maxScore = 10 * (humansMax * humansMax) * fibSum;

                if (state.score > 100000000L) System.out.println("LONG ERROR!!!! Overflow score");

                double score = (((double) state.score) / maxScore) * 10.0;
                score += 1;
                score *= score;
                score *= score;

                double humansAlive = state.humans.size() > 0 ? 10.0 : 0.0;
                humansAlive *= humansAlive;
                humansAlive *= humansAlive;
                
/*                double humansAlive = (((double)state.humans.size()) / humansMax) * 10.0;
                humansAlive+=1;
                humansAlive *=humansAlive;
                humansAlive *=humansAlive;
                
                double zombiesAlive = 10.0 - ((((double)state.zombies.size()) / zombiesMax) * 10.0);
                zombiesAlive+=1;
                zombiesAlive *=zombiesAlive;
                zombiesAlive *=zombiesAlive;*/


                double pathLength =
                        10.0 - ((((double) genes.size()) / settings.INDIVIDUAL_LENGTH) * 10.0);


                //Get fitness
                double fitnessScore = 0;
                if (state.humans.size() > 0) fitnessScore += score * settings.GENE_WEIGHTS[0];
                fitnessScore += humansAlive * settings.GENE_WEIGHTS[1];
                //if(state.humans.size()>0) fitnessScore += zombiesAlive * settings.GENE_WEIGHTS[2];
                if (state.humans.size() > 0 && state.zombies.isEmpty()) fitnessScore +=
                        pathLength * settings.GENE_WEIGHTS[3];
                return fitnessScore;
            }

            @Override
            public boolean equals(Object o) {
                if (this == o) return true;
                if (o == null || getClass() != o.getClass()) return false;

                Individual that = (Individual) o;

                if (Double.compare(that.fitness, fitness) != 0) return false;
                if (zombiesMax != that.zombiesMax) return false;
                if (humansMax != that.humansMax) return false;
                if (genes != null ? !genes.equals(that.genes) : that.genes != null) return false;
                if (states != null ? !states.equals(that.states) : that.states != null)
                    return false;
                return state != null ? state.equals(that.state) : that.state == null;
            }

            @Override
            public int hashCode() {
                int result;
                long temp;
                result = genes != null ? genes.hashCode() : 0;
                result = 31 * result + (states != null ? states.hashCode() : 0);
                result = 31 * result + (state != null ? state.hashCode() : 0);
                temp   = Double.doubleToLongBits(fitness);
                result = 31 * result + (int) (temp ^ (temp >>> 32));
                result = 31 * result + zombiesMax;
                result = 31 * result + humansMax;
                return result;
            }
        }

        public static class Gene {
            private final int[] vars;

            public Gene() {
                vars = new int[2];
                double distanceWeight = Genetic.R.nextDouble();
                //[0] = distanceWeight > 0.50 ? 1000 : (int) (2000 * distanceWeight);
                vars[0] = distanceWeight > 0.50 ? 1000 : 0;
                //vars[0] = (int) (1000 * distanceWeight);
                vars[1] = Genetic.R.nextInt(16);
            }

            public Gene(int[] vars) {
                this.vars = vars;
            }

            int[] getValues() {
                return vars;
            }

            @Override
            public boolean equals(Object o) {
                if (this == o) return true;
                if (o == null || getClass() != o.getClass()) return false;

                Gene gene = (Gene) o;

                return Arrays.equals(vars, gene.vars);
            }

            @Override
            public int hashCode() {
                return Arrays.hashCode(vars);
            }

            @Override
            public String toString() {
                return "Gene{" + "vars=" + Arrays.toString(vars) + '}';
            }
        }

        public static class Solver {
            public Settings settings;

            public Solver(Settings settings) {
                this.settings = settings;
            }

            public List<Individual> evolve(List<Individual> pop, State state, long startTime,
                                           int[] crossovers) {
                if (settings.REMOVE_STEP) removeFirstStep(pop);
                if (pop.isEmpty()) pop.addAll(
                        getRandomIndividuals(settings.DESIRED_POPULATION_SIZE));
                int counter = 0;
                while (System.currentTimeMillis() - startTime < settings.EVALUATE_TIME) {
                    int eliteSize = (int) (pop.size() * settings.ELITISM_PERCENTAGE);
                    List<Individual> elite = new ArrayList<>(pop.subList(0, eliteSize));
                    List<Individual> newPop = crossover(pop);
                    mutate(newPop);
                    elite.addAll(newPop);
                    newPop = elite;
                    newPop = resizeList(newPop);
                    for (Individual i : newPop) i.simulate(state, settings);
                    if (settings.REMOVE_DUPLICATES) newPop = new ArrayList<>(new HashSet<>(newPop));
                    pop = newPop;
                    pop.sort(Comparator.comparingDouble(Individual::getFitness).reversed());
                    counter++;
                }
                crossovers[0] += counter;

                return pop;
            }

            public void removeFirstStep(List<Individual> pop) {
                for (int i = pop.size() - 1; i >= 0; i--){
                    Individual ind = pop.get(i);
                    if (ind.getGenes().size() < 2) pop.remove(i);
                    else ind.removeFirstStep();
                }
            }

            private void mutate(List<Individual> pop) {
                for (Individual ind : pop){
                    List<Gene> genes = ind.getGenes();
                    for (int j = 0; j < genes.size(); j++)
                        if (Genetic.R.nextDouble() < settings.MUTATION_CHANCE) genes.set(j,
                                new Gene());
                }
            }

            private List<Individual> crossover(List<Individual> pop) {
                List<Individual> result = new ArrayList<>();
                while (result.size() <
                       settings.DESIRED_POPULATION_SIZE * settings.CROSSOVER_PERCENTAGE) {
                    Individual i1 = settings.selectType.doSelect(pop, settings.TOURNAMENT_SIZE);
                    Individual i2 = settings.selectType.doSelect(pop, settings.TOURNAMENT_SIZE);
                    if (settings.RANDOM_CROSSOVER_ON_DUPLICATE){
                        ArrayList<Gene> randGenes = new ArrayList<>();
                        for (int i = 0; i < settings.INDIVIDUAL_LENGTH; i++)
                             randGenes.add(new Gene());
                        i2 = new Individual(randGenes);
                    } else {
                        int counter = 0;
                        while (i1 == i2 && counter < 100) {
                            i2 = settings.selectType.doSelect(pop, settings.TOURNAMENT_SIZE);
                            counter++;
                        }
                    }
                    settings.crossoverType.doCrossover(result, i1, i2);
                }
                return result;
            }

            private List<Individual> resizeList(List<Individual> l1) {
                if (l1.size() > settings.DESIRED_POPULATION_SIZE) l1 = new ArrayList<>(
                        l1.subList(0, settings.DESIRED_POPULATION_SIZE));
                else if (l1.size() < settings.DESIRED_POPULATION_SIZE) l1.addAll(
                        getRandomIndividuals(settings.DESIRED_POPULATION_SIZE - l1.size()));
                return l1;
            }

            public List<Individual> getRandomIndividuals(int length) {
                List<Individual> result = new ArrayList<>();
                for (int i = 0; i < length; i++)
                     result.add(new Individual());
                return result;
            }
        }

        public static class Settings {
            public Genetic.SelectType selectType = SelectType.TOURNAMENT;
            public Genetic.CrossoverType crossoverType = CrossoverType.POINT;
            public boolean RANDOM_CROSSOVER_ON_DUPLICATE = true;
            public boolean REMOVE_DUPLICATES = false;
            public boolean REMOVE_STEP = false;
            public int TOURNAMENT_SIZE = 5;
            public double ELITISM_PERCENTAGE = 0.3;
            public int INDIVIDUAL_LENGTH = 50;
            public int DESIRED_POPULATION_SIZE = 50;
            public double CROSSOVER_PERCENTAGE = 0.7;
            public double MUTATION_CHANCE = 0.02;
            public double[] GENE_WEIGHTS = new double[]{1, 0.5, 0.5, 0.1};
            public long EVALUATE_TIME = 95;

            public Settings() {
            }

            public Settings(Settings s) {
                this.selectType                    = s.selectType;
                this.crossoverType                 = s.crossoverType;
                this.RANDOM_CROSSOVER_ON_DUPLICATE = s.RANDOM_CROSSOVER_ON_DUPLICATE;
                this.REMOVE_DUPLICATES             = s.REMOVE_DUPLICATES;
                this.REMOVE_STEP                   = s.REMOVE_STEP;
                this.TOURNAMENT_SIZE               = s.TOURNAMENT_SIZE;
                this.ELITISM_PERCENTAGE            = s.ELITISM_PERCENTAGE;
                this.INDIVIDUAL_LENGTH             = s.INDIVIDUAL_LENGTH;
                this.DESIRED_POPULATION_SIZE       = s.DESIRED_POPULATION_SIZE;
                this.CROSSOVER_PERCENTAGE          = s.CROSSOVER_PERCENTAGE;
                this.MUTATION_CHANCE               = s.MUTATION_CHANCE;
                this.GENE_WEIGHTS                  = s.GENE_WEIGHTS;
                this.EVALUATE_TIME                 = s.EVALUATE_TIME;
            }
        }
    }

    static class KdTree {
        private static Node insertNode(Node root, Actor actor, boolean isAlive, boolean isX) {
            if (root == null) return new Node(actor, isAlive);
            if (isX ? (actor.x < root.actor.x) : (actor.y < root.actor.y)) 
                root.left = insertNode(root.left, actor, isAlive, !isX);
            else root.right = insertNode(root.right, actor, isAlive, !isX);
            return root;
        }

        public static Node insert(Node root, Actor actor, boolean isAlive) {
            return insertNode(root, actor, isAlive, true);
        }

        public static Node nearestNeighbour(Node root, Actor actor, boolean isAlive) {
            return searchNearestNeighbour(root, actor, Integer.MAX_VALUE, root, isAlive);
        }

        public static Node searchNearestNeighbour(Node root, Actor actor, int minDist,
                                                  Node bestNode, boolean isAlive) {
            if (root == null) return bestNode;
            int distToRoot = actor.getDistanceSqrt(root.actor);
            if (distToRoot < minDist && root.isAlive == isAlive){
                minDist  = distToRoot;
                bestNode = root;
            }
            if (root.left == null)
                return searchNearestNeighbour(root.right, actor, minDist, bestNode, isAlive);
            if (root.right == null)
                return searchNearestNeighbour(root.left, actor, minDist, bestNode, isAlive);

            if (actor.getDistanceSqrt(root.left.actor) <
                actor.getDistanceSqrt(root.right.actor) && root.isAlive == isAlive) bestNode = searchNearestNeighbour(
                    root.left, actor, minDist, bestNode, isAlive);
            else bestNode = searchNearestNeighbour(root.right, actor, minDist, bestNode, isAlive);
            return bestNode;
        }

        static class Node {
            public Player.Actor actor;
            public Node left;
            public Node right;
            boolean isAlive;

            public Node(Actor actor, boolean isAlive) {
                this.isAlive = isAlive;
                this.actor   = actor;
                left         = null;
                right        = null;
            }
        }
    }
}