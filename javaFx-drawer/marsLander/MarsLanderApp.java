package marsLander;

import javafx.application.Application;
import javafx.application.Platform;
import javafx.beans.property.SimpleIntegerProperty;
import javafx.beans.property.SimpleStringProperty;
import javafx.geometry.Insets;
import javafx.scene.Node;
import javafx.scene.Scene;
import javafx.scene.canvas.Canvas;
import javafx.scene.image.Image;
import javafx.scene.image.ImageView;
import javafx.scene.layout.*;
import javafx.scene.paint.Color;
import javafx.scene.shape.Polyline;
import javafx.scene.text.Font;
import javafx.scene.text.Text;
import javafx.stage.Stage;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CountDownLatch;

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
    private List<Polyline> previousPaths = new ArrayList<>();
    
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

    public void updateObservable(int gen, int converged,int x, int y, int hSpeed, int vSpeed, int angle, int power,
                                 int fuel, Map<double[], Double> paths) {
        this.y.set(x - IMAGE_WIDTH/2);
        this.x.set(y - IMAGE_HEIGHT);
        rot.set(-angle);
        simulationText.set("generation: " + gen + " converged: " + converged);
        coordText.set("coord x: " + x + " y: " + y);
        speedText.set("speed x: " + hSpeed + " y: " + vSpeed);
        controlText.set("angle: " + angle + " power: " + power + " fuel: " + fuel);
        buildPaths(paths);
    }

    private void buildPaths(Map<double[], Double> paths) {
        Double maxScore = paths.entrySet()
                .stream()
                .filter(e->e.getValue()<10.0)
                .max(Comparator.comparingDouble(Map.Entry::getValue))
                .map(Map.Entry::getValue)
                .orElse(1.0);
        List<Polyline> newPath = new ArrayList<>();
        for (Map.Entry<double[], Double> e : paths.entrySet()){
            double[] path = e.getKey();
            Polyline polyline = new Polyline(path);
            if(e.getValue()>=10.0) polyline.setStroke(Color.WHITE);
            else if(e.getValue()<=0.0) polyline.setStroke(Color.DARKRED);
            else polyline.setStroke(Color.color(0.0, e.getValue() / maxScore ,0.0));
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
