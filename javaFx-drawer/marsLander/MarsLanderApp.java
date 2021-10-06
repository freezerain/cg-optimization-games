package marsLander;

import javafx.application.Application;
import javafx.application.Platform;
import javafx.beans.property.SimpleIntegerProperty;
import javafx.beans.property.SimpleStringProperty;
import javafx.geometry.Insets;
import javafx.scene.Node;
import javafx.scene.Scene;
import javafx.scene.canvas.Canvas;
import javafx.scene.control.Tooltip;
import javafx.scene.image.Image;
import javafx.scene.image.ImageView;
import javafx.scene.layout.*;
import javafx.scene.paint.Color;
import javafx.scene.shape.Circle;
import javafx.scene.shape.Polyline;
import javafx.scene.text.Font;
import javafx.scene.text.Text;
import javafx.stage.Stage;
import javafx.util.Duration;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.stream.DoubleStream;

public class MarsLanderApp extends Application {
    private static final CountDownLatch latch = new CountDownLatch(1);
    public static MarsLanderApp appReference = null;
    private final SimpleIntegerProperty y = new SimpleIntegerProperty();
    private final SimpleIntegerProperty x = new SimpleIntegerProperty();
    private final SimpleIntegerProperty rot = new SimpleIntegerProperty();
    private final SimpleStringProperty simulationText = new SimpleStringProperty();
    private final SimpleStringProperty coordText = new SimpleStringProperty();
    private final SimpleStringProperty speedText = new SimpleStringProperty();
    private final SimpleStringProperty controlText = new SimpleStringProperty();
    int WIDTH = 1750;
    int HEIGHT = 750;
    int IMAGE_WIDTH = 67;
    int IMAGE_HEIGHT = 62;
    private Pane pane;
    private List<Node> previousPaths = new ArrayList<>();

    public MarsLanderApp() {
        appReference = this;
    }

    public static void main(String[] args) {
        launch();
    }

    public static MarsLanderApp waitForStart() {
        try {
            latch.await();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        return appReference;
    }

    public void launchApp() {
        launch();
    }

    @Override
    public void start(Stage stage) {
        stage.setTitle("Mars Lander");
        pane = new Pane();
        pane.setBackground(new Background(
                new BackgroundFill(Color.color(0.0, 0.0, 0.0), new CornerRadii(0), new Insets(0))));
        Canvas canvas = new CanvasGridGenerator(WIDTH, HEIGHT, 10).generate(0.25,
                Color.color(1, 1, 1, 0.01));
        pane.getChildren().add(canvas);

        Text text0 = new Text();
        text0.setFill(Color.WHITE);
        text0.setFont(Font.font(18.0));
        text0.textProperty().bind(simulationText);
        Text text1 = new Text();
        text1.setFill(Color.WHITE);
        text1.setFont(Font.font(18.0));
        text1.textProperty().bind(coordText);
        Text text2 = new Text();
        text2.setFill(Color.WHITE);
        text2.setFont(Font.font(18.0));
        text2.textProperty().bind(speedText);
        Text text3 = new Text();
        text3.setFill(Color.WHITE);
        text3.setFont(Font.font(18.0));
        text3.textProperty().bind(controlText);
        VBox vbox = new VBox(text0, text1, text2, text3);
        vbox.setLayoutX(100);
        vbox.setLayoutY(100);
        pane.getChildren().add(vbox);

        Image image = new Image(getClass().getResource("MarsLander1.png").toString());
        ImageView im = new ImageView(image);
        im.layoutXProperty().bind(y);
        im.layoutYProperty().bind(x);
        im.rotateProperty().bind(rot);
        pane.getChildren().add(im);
        stage.setWidth(WIDTH);
        stage.setHeight(HEIGHT + 20);
        stage.setScene(new Scene(pane));
        stage.show();
        latch.countDown();
    }

    public void updateObservable(int gen, List<LanderPath> paths) {
        if(paths.isEmpty()) return;
        Lander lastState = paths.get(0).firstState;
        this.y.set(((int) (Math.round(lastState.x / 4.0))) - IMAGE_WIDTH / 2);
        this.x.set(((int) (Math.round(750 - lastState.y / 4.0))) - IMAGE_HEIGHT);
        rot.set(-lastState.angle);
        long safeLandedCount = paths.stream().filter(p -> p.isSafeLanded).count();
        simulationText.set("generation: " + gen + " converged: " + safeLandedCount);
        coordText.set("coord x: " + Math.round(lastState.x / 4.0) + " y: " +
                      Math.round(lastState.y / 4.0));
        speedText.set("speed x: " + ((int) Math.round(lastState.hSpeed)) + " y: " +
                      ((int) Math.round(lastState.vSpeed)));
        controlText.set("angle: " + lastState.angle + " power: " + lastState.power + " fuel: " +
                        lastState.fuel);
        buildPaths(paths);
    }

    private void buildPaths(List<LanderPath> landerPath) {
        double maxScore = landerPath.stream()
                .filter(l -> l.isLanded && !l.isSafeLanded)
                .mapToDouble(l -> l.score)
                .max()
                .orElse(1.0);
        double minScore = landerPath.stream()
                .filter(l -> l.isLanded && !l.isSafeLanded)
                .mapToDouble(l -> l.score)
                .min()
                .orElse(0.0);
        List<Node> newPath = new ArrayList<>();
        for (LanderPath lp : landerPath){
            double[] path = lp.landerList.stream()
                    .flatMapToDouble(l -> DoubleStream.of(l.x / 4.0, 750.0 - l.y / 4.0))
                    .toArray();
            Polyline polyline = new Polyline(path);
            if (lp.isSafeLanded) polyline.setStroke(Color.WHITE);
            else if (!lp.isLanded) polyline.setStroke(Color.DARKRED);
            else if (lp.score < 0) polyline.setStroke(Color.DIMGREY);
            else polyline.setStroke(
                        Color.color(0.0,  (lp.score - minScore) / (maxScore - minScore), 0.0));


            Circle c = new Circle();
            c.setRadius(5);
            c.setFill(Color.color(0.35, 0.35, 0.35));
            c.setLayoutX(lp.lastState.x / 4.0);
            c.setLayoutY(750.0 - lp.lastState.y / 4.0);
            String tp = "x: " + lp.lastState.x + " y: " + lp.lastState.y + "\n" + "hS: " +
                        lp.lastState.hSpeed + " vS: " + lp.lastState.vSpeed + "\n" + "angle: " +
                        lp.lastState.angle + " power: " + lp.lastState.power + " fuel: " +
                        lp.lastState.fuel + "\n" + "isLanded: " + lp.lastState.isLanded +
                        " isSafe: " + lp.lastState.isSafeLanded + "\nfitness: " + lp.score;
            Tooltip tooltip = new Tooltip(tp);
            tooltip.setShowDelay(Duration.millis(10));
            tooltip.setShowDuration(Duration.seconds(30));
            Tooltip.install(c, tooltip);
            newPath.add(c);

            newPath.add(polyline);
        }
        System.out.println("polyline size:" + newPath.size());
        Platform.runLater(() -> {
            pane.getChildren().removeAll(previousPaths);
            pane.getChildren().addAll(newPath);
            previousPaths = newPath;
        });
    }

    public void buildLandscape(double[] landscape) {
        Polyline land = new Polyline(landscape);
        land.setStrokeWidth(1);
        land.setStroke(Color.DARKRED);
        Platform.runLater(() -> pane.getChildren().add(land));
    }
}
