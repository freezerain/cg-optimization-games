package CodeVsZombies;

import javafx.application.Application;
import javafx.application.Platform;
import javafx.geometry.Insets;
import javafx.scene.Scene;
import javafx.scene.layout.Background;
import javafx.scene.layout.BackgroundFill;
import javafx.scene.layout.CornerRadii;
import javafx.scene.layout.Pane;
import javafx.scene.paint.Color;
import javafx.scene.shape.Circle;
import javafx.scene.shape.Polyline;
import javafx.stage.Stage;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CountDownLatch;

//16000 / 10 = 1600
//9000 / 10 = 900

public class CodeVsZombiesApp extends Application {
    private static final CountDownLatch latch = new CountDownLatch(1);
    public static CodeVsZombiesApp appReference = null;
    private Pane pane;

    public CodeVsZombiesApp() {
        appReference = this;
    }

    public static void main(String[] args) {
        launch();
    }

    public static CodeVsZombiesApp waitForStart() {
        try {
            latch.await();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        return appReference;
    }

    public void updateScene(int[] player, int[] humans, int[] zombies, int[][] paths,
                            double[] scores) {
        List<Polyline> polylines = buildLines(paths, scores);
        List<Circle> p = buildCircles(player, Color.WHITE);
        List<Circle> h = buildCircles(humans, Color.BLUE);
        List<Circle> z = buildCircles(zombies, Color.DARKORANGE);
        Platform.runLater(() -> {
            pane.getChildren().clear();
            pane.getChildren().addAll(polylines);
            pane.getChildren().addAll(p);
            pane.getChildren().addAll(h);
            pane.getChildren().addAll(z);
        });
    }

    private List<Circle> buildCircles(int[] coord, Color color) {
        var list = new ArrayList<Circle>();
        for (int i = 0; i < coord.length; i += 2){
            Circle circle = new Circle();
            circle.setRadius(25);
            circle.setFill(color);
            circle.setLayoutX(coord[i] / 10.0);
            circle.setLayoutY((coord[(i + 1)] / 10.0));
            list.add(circle);
        }
        return list;
    }

    private List<Polyline> buildLines(int[][] paths, double[] scores) {
        double minScore = Double.MAX_VALUE;
        double maxScore = Double.MIN_VALUE;
        for (double score : scores){
            minScore = Math.min(minScore, score);
            maxScore = Math.max(maxScore, score);
        }
        var list = new ArrayList<Polyline>();
        for (int i = 0; i < paths.length; i++){
            int[] path = paths[i];
            var doublePath = new double[path.length];
            for (int j = 0; j < path.length; j++){
                doublePath[j] = path[j] / 10.0;
            }
            Polyline pl = new Polyline(doublePath);
            Color color;
            if (scores[i] < 0) color = Color.DARKRED;
            else color = Color.color(0.0, (scores[i] - minScore) / (maxScore - minScore), 0.0);
            pl.setStroke(color);
            list.add(pl);
        }
        return list;
    }

    @Override
    public void start(Stage stage) throws Exception {
        stage.setTitle("Mars Lander");
        pane = new Pane();
        pane.setBackground(new Background(
                new BackgroundFill(Color.color(0.0, 0.0, 0.0), new CornerRadii(0), new Insets(0))));
        stage.setWidth(1600);
        stage.setHeight(900);
        stage.setScene(new Scene(pane));
        stage.show();
        latch.countDown();
    }


}
