package src;

import javafx.application.Application;
import javafx.scene.Scene;
import javafx.scene.canvas.Canvas;
import javafx.scene.image.Image;
import javafx.scene.layout.*;
import javafx.scene.paint.Color;
import javafx.scene.shape.Circle;
import javafx.scene.shape.Line;
import javafx.stage.Stage;

public class Mars_Lander_App extends Application {
    public static void main(String[] args) {
        launch();
    }

    @Override
    public void start(Stage stage) throws Exception {
        stage.setTitle("Mars Lander");
        //stage.initStyle(StageStyle.UNIFIED);
        Pane pane = new Pane();
        //BackgroundFill b1 = new BackgroundFill(Color.color(0.0, 0.0, 0.0), new CornerRadii(0),
              //  new Insets(0));
        //BackgroundFill b2 = new BackgroundFill(Color.color(0.5, 0.5, 0.5), new CornerRadii(5),
            //    new Insets(0));
       // pane.setBackground(new Background(b1, b2));
        Image background = new Image("res/MarsLander1.jpg");
        pane.setBackground(new Background(
                new BackgroundImage(background, BackgroundRepeat.REPEAT, BackgroundRepeat.NO_REPEAT,
                        BackgroundPosition.DEFAULT, BackgroundSize.DEFAULT)));
        Canvas canvas = new CanvasGridGenerator(1800, 1000, 10).generate(0.25, Color.color(1, 1, 1, 0.01));
        Circle circle = new Circle();
        circle.setFill(Color.color(1,1,1));
        circle.setRadius(20);
        circle.setLayoutX(700);
        circle.setLayoutY(200);
        pane.getChildren().add(circle);
        Line line = new Line(500, 200, 750, 600);
        line.setStrokeWidth(5);
        pane.getChildren().add(line);
        Line line2 = new Line(450, 150, 700, 550);
        line2.setStrokeWidth(10);
        pane.getChildren().add(line2);
        pane.getChildren().add(canvas);
        stage.setWidth(1800);
        stage.setHeight(1000);
        stage.setScene(new Scene(pane));
        stage.show();
    }
}
