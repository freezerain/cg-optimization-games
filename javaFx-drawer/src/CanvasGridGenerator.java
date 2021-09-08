package src;

import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.paint.Color;
import javafx.scene.paint.Paint;

public class CanvasGridGenerator {
    double width;
    double height;
    double step;

    public CanvasGridGenerator(double width, double height, double step) {
        this.width  = width;
        this.height = height;
        this.step   = step;
    }

    public Canvas generate(double lineWidth, Paint color){
        Canvas canvas = new Canvas();
        canvas.prefWidth(Double.MAX_VALUE);
        canvas.prefHeight(Double.MAX_VALUE);
        canvas.setWidth(width);
        canvas.setHeight(height);
        GraphicsContext gc = canvas.getGraphicsContext2D();
        gc.setLineWidth(lineWidth);
        gc.setStroke(color);
        for (double i = 0.5; i < width; i+=step){
            gc.moveTo(i, 0);
            gc.lineTo(i, height);
            gc.stroke();
        }
        for (double i = 0.5; i < height; i+=step){
            gc.moveTo(0, i);
            gc.lineTo(width, i);
            gc.stroke();
        }
        return canvas;
    }
}
