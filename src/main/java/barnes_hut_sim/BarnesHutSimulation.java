package barnes_hut_sim;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * BarnesHutSimulation – A complex n-body simulation using the Barnes-Hut algorithm
 * to approximate gravitational forces efficiently. A quadtree is built so that
 * for each body's force calculation, either the exact force is used or an approximation
 * if the body is sufficiently far away.
 */
public class BarnesHutSimulation extends JPanel implements ActionListener, MouseListener, MouseMotionListener {

    // --- Simulation Constants ---
    private static final int WIDTH = 800;
    private static final int HEIGHT = 800;
    private static final double DT = 0.1;        // Time step
    private static final double G = 1.0;         // Scaled gravitational constant
    private static final double THETA = 0.5;     // Opening angle for Barnes-Hut (smaller => more accurate)
    private static final int NUM_BODIES = 500;   // Initial number of bodies

    private final transient List<Body> bodies;

    // For mouse interactions (a mouse click inserts a new body)
    public BarnesHutSimulation() {
        setPreferredSize(new Dimension(WIDTH, HEIGHT));
        setBackground(Color.BLACK);
        bodies = new ArrayList<>();
        Random rand = new Random();
        for (int i = 0; i < NUM_BODIES; i++) {
            double x = rand.nextDouble() * WIDTH;
            double y = rand.nextDouble() * HEIGHT;
            double vx = (rand.nextDouble() - 0.5) * 2;
            double vy = (rand.nextDouble() - 0.5) * 2;
            double mass = rand.nextDouble() * 10 + 1;
            bodies.add(new Body(x, y, vx, vy, mass));
        }
        addMouseListener(this);
        addMouseMotionListener(this);
        Timer timer = new Timer(16, this); // ~60 FPS
        timer.start();
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        updateSimulation();
        repaint();
    }

    /**
     * Updates the simulation:
     * - Defines a Quad that covers the entire simulation area.
     * - Creates a quadtree (Node) and inserts all bodies into it.
     * - For each body, uses the tree to compute the resulting force.
     * - Finally, updates each body's velocity and position.
     */
    private void updateSimulation() {
        Quad quad = new Quad(0, 0, WIDTH);
        Node tree = new Node(quad);
        for (Body b : bodies) {
            tree.insert(b);
        }
        for (Body b : bodies) {
            b.resetForce();
            tree.updateForce(b);
        }
        for (Body b : bodies) {
            b.update(DT);
        }
    }

    @Override
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);
        Graphics2D g2d = (Graphics2D) g;
        // Draw all bodies as small circles in white
        g2d.setColor(Color.WHITE);
        for (Body b : bodies) {
            int r = (int) Math.max(2, Math.cbrt(b.mass)); // Radius proportional to the cube root of mass
            g2d.fillOval((int) (b.x - r), (int) (b.y - r), 2 * r, 2 * r);
        }
    }

    // --- Mouse interaction: a mouse click adds a new body ---
    @Override
    public void mousePressed(MouseEvent e) {
        double x = e.getX();
        double y = e.getY();
        // Add a body with zero velocity and a fixed mass
        bodies.add(new Body(x, y, 0, 0, 20));
    }

    @Override
    public void mouseReleased(MouseEvent e) { /* Not used */ }

    @Override
    public void mouseDragged(MouseEvent e) { /* Not used */ }

    @Override
    public void mouseMoved(MouseEvent e) { /* Not used */ }

    @Override
    public void mouseClicked(MouseEvent e) { /* Not used */ }

    @Override
    public void mouseEntered(MouseEvent e) { /* Not used */ }

    @Override
    public void mouseExited(MouseEvent e) { /* Not used */ }

    // --------------------- Inner Classes ---------------------

    /**
     * Body – Represents a body in the simulation with position, velocity, force, and mass.
     */
    private static class Body {
        double x;    // Position x
        double y;    // Position y
        double vx;    // Velocity x
        double vy;    // Velocity y
        double fx;    // Accumulated force x
        double fy;    // Accumulated force y
        double mass;    // Mass

        Body(double x, double y, double vx, double vy, double mass) {
            this.x = x;
            this.y = y;
            this.vx = vx;
            this.vy = vy;
            this.mass = mass;
        }

        // Resets the accumulated forces
        public void resetForce() {
            fx = 0;
            fy = 0;
        }

        /**
         * Adds the force exerted by body b on this body.
         * Uses a softening parameter to avoid singularities.
         */
        public void addForce(Body b) {
            double dx = b.x - this.x;
            double dy = b.y - this.y;
            double dist = Math.sqrt(dx * dx + dy * dy);
            double eps = 3; // Softening factor
            double force = (G * this.mass * b.mass) / ((dist * dist) + (eps * eps));
            // Normalize direction vectors and multiply by force
            fx += force * dx / dist;
            fy += force * dy / dist;
        }

        // Updates velocity and position based on accumulated force
        public void update(double dt) {
            vx += (fx / mass) * dt;
            vy += (fy / mass) * dt;
            x += vx * dt;
            y += vy * dt;
            // Reflect from the borders
            if (x < 0) {
                x = 0;
                vx = -vx;
            }
            if (x > WIDTH) {
                x = WIDTH;
                vx = -vx;
            }
            if (y < 0) {
                y = 0;
                vy = -vy;
            }
            if (y > HEIGHT) {
                y = HEIGHT;
                vy = -vy;
            }
        }
    }

    /**
     * Quad – Represents a square region in the simulation.
     */
    private static class Quad {
        double x;
        double y;
        double length;     // Side length of the square

        Quad(double x, double y, double length) {
            this.x = x;
            this.y = y;
            this.length = length;
        }

        // Checks if a point (cx, cy) lies inside this square
        public boolean contains(double cx, double cy) {
            return (cx >= x && cx <= x + length && cy >= y && cy <= y + length);
        }

        public Quad nw() {
            return new Quad(x, y, length / 2);
        }

        public Quad ne() {
            return new Quad(x + length / 2, y, length / 2);
        }

        public Quad sw() {
            return new Quad(x, y + length / 2, length / 2);
        }

        public Quad se() {
            return new Quad(x + length / 2, y + length / 2, length / 2);
        }
    }

    /**
     * Node – Represents a node in the quadtree.
     * A node can be either external (containing a single body)
     * or internal (containing multiple bodies and subdivided into quadrants).
     */
    private static class Node {
        Quad quad;             // The region this node covers
        Body body;             // Holds a body if this node is external
        Node nw;               // Node
        Node ne;               // Node
        Node sw;               // Node
        Node se;               // Node
        double mass;           // Total mass in this node
        double centerX; // Center of mass
        double centerY; // Center of mass

        Node(Quad quad) {
            this.quad = quad;
            this.body = null;
            this.mass = 0;
            this.centerX = 0;
            this.centerY = 0;
        }

        // Checks if this node is external (a leaf in the quadtree)
        public boolean isExternal() {
            return (nw == null && ne == null && sw == null && se == null);
        }

        /**
         * Inserts a body into this node (and its child nodes), updating
         * the total mass and center of mass.
         */
        public void insert(Body b) {
            if (!quad.contains(b.x, b.y)) {
                return; // Body does not lie in this region
            }
            if (body == null && isExternal()) {
                // Empty external node: place the body here
                body = b;
                mass = b.mass;
                centerX = b.x;
                centerY = b.y;
            } else {
                // If the node is external, subdivide and move the existing body
                if (isExternal()) {
                    subdivide();
                    if (body != null) {
                        putBody(body);
                        body = null;
                    }
                }
                // Insert the new body into the correct child node
                putBody(b);
                // Update total mass and center of mass
                double totalMass = mass + b.mass;
                centerX = (centerX * mass + b.x * b.mass) / totalMass;
                centerY = (centerY * mass + b.y * b.mass) / totalMass;
                mass = totalMass;
            }
        }

        // Subdivides this node into four equally sized child nodes
        private void subdivide() {
            nw = new Node(quad.nw());
            ne = new Node(quad.ne());
            sw = new Node(quad.sw());
            se = new Node(quad.se());
        }

        // Places a body into the appropriate child node
        private void putBody(Body b) {
            if (quad.nw().contains(b.x, b.y)) {
                nw.insert(b);
            } else if (quad.ne().contains(b.x, b.y)) {
                ne.insert(b);
            } else if (quad.sw().contains(b.x, b.y)) {
                sw.insert(b);
            } else if (quad.se().contains(b.x, b.y)) {
                se.insert(b);
            }
        }

        /**
         * Updates the force on a body b by traversing the quadtree
         * and using either the entire node (if sufficiently far away)
         * or descending to child nodes for a more detailed calculation.
         */
        public void updateForce(Body b) {
            if (body == null && mass == 0) return; // No contribution
            if (isExternal()) {
                if (body != b) {
                    b.addForce(body);
                }
            } else {
                double dx = centerX - b.x;
                double dy = centerY - b.y;
                double dist = Math.sqrt(dx * dx + dy * dy);
                if (quad.length / dist < THETA) {
                    // Approximation: treat this node as a single body
                    Body pseudo = new Body(centerX, centerY, 0, 0, mass);
                    b.addForce(pseudo);
                } else {
                    if (nw != null) nw.updateForce(b);
                    if (ne != null) ne.updateForce(b);
                    if (sw != null) sw.updateForce(b);
                    if (se != null) se.updateForce(b);
                }
            }
        }
    }

    // --------------------- Main Method ---------------------
    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            JFrame frame = new JFrame("Barnes-Hut Simulation");
            frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
            BarnesHutSimulation sim = new BarnesHutSimulation();
            frame.add(sim);
            frame.pack();
            frame.setLocationRelativeTo(null);
            frame.setVisible(true);
        });
    }
}
