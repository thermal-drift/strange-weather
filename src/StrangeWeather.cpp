/*
 * Strange Weather - Chaotic CV Generator
 * 
 * Three independent strange attractor banks producing 16 CV outputs.
 * Based on Lorenz, Rössler, Thomas, and Halvorsen attractors.
 */

#include "plugin.hpp"
#include <cmath>

// Attractor types
enum AttractorType {
    LORENZ = 0,
    ROSSLER = 1,
    THOMAS = 2,
    HALVORSEN = 3
};

// Attractor state and parameters
struct Attractor {
    double x, y, z;
    AttractorType type;
    float chaos = 0.5f; // 0-1, controls primary chaos parameter

    // Bounding box tracking for normalization
    double minX, maxX, minY, maxY, minZ, maxZ;

    Attractor() {
        // Initialize with slightly random starting point
        x = 0.1 + (random::uniform() - 0.5) * 0.1;
        y = 0.0 + (random::uniform() - 0.5) * 0.1;
        z = 0.0 + (random::uniform() - 0.5) * 0.1;
        type = LORENZ;

        // Initial bounds (will adapt)
        minX = -20.0; maxX = 20.0;
        minY = -30.0; maxY = 30.0;
        minZ = 0.0;   maxZ = 50.0;
    }

    // Compute derivatives for current state (chaos affects primary parameter)
    void derivatives(double& dx, double& dy, double& dz) {
        switch (type) {
            case LORENZ: {
                // σ = 10, β = 8/3, ρ varies with chaos: 20-30
                const double sigma = 10.0;
                const double rho = 20.0 + chaos * 10.0;  // chaos: periodic → chaotic
                const double beta = 8.0 / 3.0;
                dx = sigma * (y - x);
                dy = x * (rho - z) - y;
                dz = x * y - beta * z;
                break;
            }
            case ROSSLER: {
                // a = 0.2, b = 0.2, c varies with chaos: 4-7
                const double a = 0.2;
                const double b = 0.2;
                const double c = 4.0 + chaos * 3.0;  // chaos: tight spiral → wild
                dx = -y - z;
                dy = x + a * y;
                dz = b + z * (x - c);
                break;
            }
            case THOMAS: {
                // b varies with chaos: 0.3-0.15 (inverted - lower = more chaos)
                const double b = 0.3 - chaos * 0.15;  // chaos: damped → sustained
                dx = std::sin(y) - b * x;
                dy = std::sin(z) - b * y;
                dz = std::sin(x) - b * z;
                break;
            }
            case HALVORSEN: {
                // a varies with chaos: 1.4-2.0
                const double a = 1.4 + chaos * 0.6;  // chaos: mild → aggressive
                dx = -a * x - 4.0 * y - 4.0 * z - y * y;
                dy = -a * y - 4.0 * z - 4.0 * x - z * z;
                dz = -a * z - 4.0 * x - 4.0 * y - x * x;
                break;
            }
        }
    }
    
    // RK4 integration step
    void step(double dt) {
        double k1x, k1y, k1z;
        double k2x, k2y, k2z;
        double k3x, k3y, k3z;
        double k4x, k4y, k4z;
        
        double ox = x, oy = y, oz = z;
        
        // k1
        derivatives(k1x, k1y, k1z);
        
        // k2
        x = ox + 0.5 * dt * k1x;
        y = oy + 0.5 * dt * k1y;
        z = oz + 0.5 * dt * k1z;
        derivatives(k2x, k2y, k2z);
        
        // k3
        x = ox + 0.5 * dt * k2x;
        y = oy + 0.5 * dt * k2y;
        z = oz + 0.5 * dt * k2z;
        derivatives(k3x, k3y, k3z);
        
        // k4
        x = ox + dt * k3x;
        y = oy + dt * k3y;
        z = oz + dt * k3z;
        derivatives(k4x, k4y, k4z);
        
        // Final update
        x = ox + (dt / 6.0) * (k1x + 2.0 * k2x + 2.0 * k3x + k4x);
        y = oy + (dt / 6.0) * (k1y + 2.0 * k2y + 2.0 * k3y + k4y);
        z = oz + (dt / 6.0) * (k1z + 2.0 * k2z + 2.0 * k3z + k4z);
        
        // Update bounding box with slow decay toward current values
        const double decay = 0.9999;
        minX = std::min(minX * decay + x * (1.0 - decay), x);
        maxX = std::max(maxX * decay + x * (1.0 - decay), x);
        minY = std::min(minY * decay + y * (1.0 - decay), y);
        maxY = std::max(maxY * decay + y * (1.0 - decay), y);
        minZ = std::min(minZ * decay + z * (1.0 - decay), z);
        maxZ = std::max(maxZ * decay + z * (1.0 - decay), z);
    }
    
    // Get normalized outputs (-5V to +5V)
    float getNormX() {
        double range = maxX - minX;
        if (range < 0.001) range = 0.001;
        return (float)(((x - minX) / range) * 10.0 - 5.0);
    }
    
    float getNormY() {
        double range = maxY - minY;
        if (range < 0.001) range = 0.001;
        return (float)(((y - minY) / range) * 10.0 - 5.0);
    }
    
    float getNormZ() {
        double range = maxZ - minZ;
        if (range < 0.001) range = 0.001;
        return (float)(((z - minZ) / range) * 10.0 - 5.0);
    }
};


struct StrangeWeather : Module {
    enum ParamIds {
        RATE_A_PARAM,
        RATE_B_PARAM,
        RATE_C_PARAM,
        RATE_D_PARAM,
        SHAPE_A_PARAM,
        SHAPE_B_PARAM,
        SHAPE_C_PARAM,
        SHAPE_D_PARAM,
        RANGE_A_PARAM,
        RANGE_B_PARAM,
        RANGE_C_PARAM,
        RANGE_D_PARAM,
        VOLTAGE_A_PARAM,
        VOLTAGE_B_PARAM,
        VOLTAGE_C_PARAM,
        VOLTAGE_D_PARAM,
        CHAOS_A_PARAM,
        CHAOS_B_PARAM,
        CHAOS_C_PARAM,
        CHAOS_D_PARAM,
        TRAIL_PARAM,
        NUM_PARAMS
    };
    
    enum InputIds {
        NUM_INPUTS
    };
    
    enum OutputIds {
        // Bank A
        A_X_OUTPUT,
        A_Y_OUTPUT,
        A_Z_OUTPUT,
        A_SUM_OUTPUT,
        // Bank B
        B_X_OUTPUT,
        B_Y_OUTPUT,
        B_Z_OUTPUT,
        B_SUM_OUTPUT,
        // Bank C
        C_X_OUTPUT,
        C_Y_OUTPUT,
        C_Z_OUTPUT,
        C_SUM_OUTPUT,
        // Bank D
        D_X_OUTPUT,
        D_Y_OUTPUT,
        D_Z_OUTPUT,
        D_SUM_OUTPUT,
        // Combined
        COMB_SUM_OUTPUT,
        COMB_RECT_OUTPUT,
        COMB_INV_OUTPUT,
        COMB_DIST_OUTPUT,
        NUM_OUTPUTS
    };
    
    enum LightIds {
        NUM_LIGHTS
    };
    
    // Four attractor banks
    Attractor attractors[4];

    // Smoothed outputs (one-pole lowpass)
    float smoothedX[4] = {0.f, 0.f, 0.f, 0.f};
    float smoothedY[4] = {0.f, 0.f, 0.f, 0.f};
    float smoothedZ[4] = {0.f, 0.f, 0.f, 0.f};

    // Display state
    int displayMode = 5; // 0=A, 1=B, 2=C, 3=D, 4=Combined, 5=All, 6=Ajman (if enabled)
    bool display3D = false; // Toggle between 2D and 3D view
    int displayStyle = 0; // 0=Trace, 1=Lissajous (tiny dots), 2=Scope
    bool ajmanEnabled = false; // Easter egg: enables Ajman mode in cycle

    // Trail history for display (ring buffer)
    static const int MAX_TRAIL_LENGTH = 4096;
    float trailX[4][MAX_TRAIL_LENGTH] = {};
    float trailY[4][MAX_TRAIL_LENGTH] = {};
    float trailZ[4][MAX_TRAIL_LENGTH] = {};
    int trailIndex = 0;

    // Combined trail
    float combTrailX[MAX_TRAIL_LENGTH] = {};
    float combTrailY[MAX_TRAIL_LENGTH] = {};
    float combTrailZ[MAX_TRAIL_LENGTH] = {};

    // Sample counter for trail updates
    int trailCounter = 0;
    
    StrangeWeather() {
        config(NUM_PARAMS, NUM_INPUTS, NUM_OUTPUTS, NUM_LIGHTS);

        // Rate knobs (fine control within selected range)
        configParam(RATE_A_PARAM, 0.f, 1.f, 0.5f, "Rate A");
        configParam(RATE_B_PARAM, 0.f, 1.f, 0.5f, "Rate B");
        configParam(RATE_C_PARAM, 0.f, 1.f, 0.5f, "Rate C");
        configParam(RATE_D_PARAM, 0.f, 1.f, 0.5f, "Rate D");

        // Range switches (0=Low, 1=Med, 2=High)
        configSwitch(RANGE_A_PARAM, 0.f, 2.f, 1.f, "Range A", {"Low (5-20 min)", "Med (1s-2min)", "High (0.1-10s)"});
        configSwitch(RANGE_B_PARAM, 0.f, 2.f, 1.f, "Range B", {"Low (5-20 min)", "Med (1s-2min)", "High (0.1-10s)"});
        configSwitch(RANGE_C_PARAM, 0.f, 2.f, 1.f, "Range C", {"Low (5-20 min)", "Med (1s-2min)", "High (0.1-10s)"});
        configSwitch(RANGE_D_PARAM, 0.f, 2.f, 1.f, "Range D", {"Low (5-20 min)", "Med (1s-2min)", "High (0.1-10s)"});

        // Shape switches (0-3)
        configSwitch(SHAPE_A_PARAM, 0.f, 3.f, 0.f, "Shape A", {"Lorenz", "Rössler", "Thomas", "Halvorsen"});
        configSwitch(SHAPE_B_PARAM, 0.f, 3.f, 1.f, "Shape B", {"Lorenz", "Rössler", "Thomas", "Halvorsen"});
        configSwitch(SHAPE_C_PARAM, 0.f, 3.f, 2.f, "Shape C", {"Lorenz", "Rössler", "Thomas", "Halvorsen"});
        configSwitch(SHAPE_D_PARAM, 0.f, 3.f, 3.f, "Shape D", {"Lorenz", "Rössler", "Thomas", "Halvorsen"});

        // Voltage switches (0=±5V, 1=±10V, 2=0-5V, 3=0-10V)
        configSwitch(VOLTAGE_A_PARAM, 0.f, 3.f, 0.f, "Voltage A", {"±5V", "±10V", "0-5V", "0-10V"});
        configSwitch(VOLTAGE_B_PARAM, 0.f, 3.f, 0.f, "Voltage B", {"±5V", "±10V", "0-5V", "0-10V"});
        configSwitch(VOLTAGE_C_PARAM, 0.f, 3.f, 0.f, "Voltage C", {"±5V", "±10V", "0-5V", "0-10V"});
        configSwitch(VOLTAGE_D_PARAM, 0.f, 3.f, 0.f, "Voltage D", {"±5V", "±10V", "0-5V", "0-10V"});

        // Chaos knobs (0-1)
        configParam(CHAOS_A_PARAM, 0.f, 1.f, 0.5f, "Chaos A", "%", 0.f, 100.f);
        configParam(CHAOS_B_PARAM, 0.f, 1.f, 0.5f, "Chaos B", "%", 0.f, 100.f);
        configParam(CHAOS_C_PARAM, 0.f, 1.f, 0.5f, "Chaos C", "%", 0.f, 100.f);
        configParam(CHAOS_D_PARAM, 0.f, 1.f, 0.5f, "Chaos D", "%", 0.f, 100.f);

        // Trail length knob (64-4096)
        configParam(TRAIL_PARAM, 0.f, 1.f, 0.7f, "Trail Length", "", 0.f, 1.f);

        // Output labels
        configOutput(A_X_OUTPUT, "Bank A X");
        configOutput(A_Y_OUTPUT, "Bank A Y");
        configOutput(A_Z_OUTPUT, "Bank A Z");
        configOutput(A_SUM_OUTPUT, "Bank A Sum");
        configOutput(B_X_OUTPUT, "Bank B X");
        configOutput(B_Y_OUTPUT, "Bank B Y");
        configOutput(B_Z_OUTPUT, "Bank B Z");
        configOutput(B_SUM_OUTPUT, "Bank B Sum");
        configOutput(C_X_OUTPUT, "Bank C X");
        configOutput(C_Y_OUTPUT, "Bank C Y");
        configOutput(C_Z_OUTPUT, "Bank C Z");
        configOutput(C_SUM_OUTPUT, "Bank C Sum");
        configOutput(D_X_OUTPUT, "Bank D X");
        configOutput(D_Y_OUTPUT, "Bank D Y");
        configOutput(D_Z_OUTPUT, "Bank D Z");
        configOutput(D_SUM_OUTPUT, "Bank D Sum");
        configOutput(COMB_SUM_OUTPUT, "Combined Sum");
        configOutput(COMB_RECT_OUTPUT, "Combined Rectified");
        configOutput(COMB_INV_OUTPUT, "Combined Inverted");
        configOutput(COMB_DIST_OUTPUT, "Combined Inverse Distance");
    }
    
    void cycleDisplay() {
        int maxModes = ajmanEnabled ? 7 : 6;
        displayMode = (displayMode + 1) % maxModes;
    }

    void cycleDisplayStyle() {
        displayStyle = (displayStyle + 1) % 3;
    }

    // Get current trail length from param (64 to 4096)
    int getTrailLength() {
        float knob = params[TRAIL_PARAM].getValue();
        // Exponential scaling: 64 at 0, 4096 at 1
        return (int)(64.f * std::pow(64.f, knob));
    }

    // Helper to calculate rate from range and knob
    float calculateRate(int range, float knob) {
        float minHz, maxHz;
        switch (range) {
            case 0: minHz = 0.0008f; maxHz = 0.003f; break;  // Low: 5-20 min
            case 1: minHz = 0.008f;  maxHz = 1.0f;   break;  // Med: 1s-2min
            default: minHz = 0.1f;   maxHz = 10.0f;  break;  // High: 0.1-10s
        }
        return minHz * std::pow(maxHz / minHz, knob);
    }

    // Helper to scale output based on voltage mode
    float scaleVoltage(float normalized, int voltageMode) {
        // normalized is -1 to +1
        switch (voltageMode) {
            case 0: return normalized * 5.0f;                    // ±5V
            case 1: return normalized * 10.0f;                   // ±10V
            case 2: return (normalized + 1.0f) * 2.5f;           // 0-5V
            case 3: return (normalized + 1.0f) * 5.0f;           // 0-10V
            default: return normalized * 5.0f;
        }
    }

    void process(const ProcessArgs& args) override {
        // Smoothing coefficient (lower = smoother, ~0.001 at 48kHz gives nice smooth output)
        float smoothCoeff = 6.0f / args.sampleRate;  // ~125ms time constant

        // Process each bank
        int rangeParams[4] = {RANGE_A_PARAM, RANGE_B_PARAM, RANGE_C_PARAM, RANGE_D_PARAM};
        int rateParams[4] = {RATE_A_PARAM, RATE_B_PARAM, RATE_C_PARAM, RATE_D_PARAM};
        int shapeParams[4] = {SHAPE_A_PARAM, SHAPE_B_PARAM, SHAPE_C_PARAM, SHAPE_D_PARAM};
        int voltageParams[4] = {VOLTAGE_A_PARAM, VOLTAGE_B_PARAM, VOLTAGE_C_PARAM, VOLTAGE_D_PARAM};
        int chaosParams[4] = {CHAOS_A_PARAM, CHAOS_B_PARAM, CHAOS_C_PARAM, CHAOS_D_PARAM};

        float bankOutputs[4][4]; // [bank][x,y,z,sum]

        for (int i = 0; i < 4; i++) {
            // Get parameters
            int range = (int)params[rangeParams[i]].getValue();
            float rateKnob = params[rateParams[i]].getValue();
            int voltageMode = (int)params[voltageParams[i]].getValue();

            // Set attractor type and chaos
            attractors[i].type = (AttractorType)(int)params[shapeParams[i]].getValue();
            attractors[i].chaos = params[chaosParams[i]].getValue();

            // Calculate rate
            float rate = calculateRate(range, rateKnob);

            // Adaptive time step
            float dt = rate / args.sampleRate;
            const float maxDt = 0.01f;
            int steps = (int)std::ceil(dt / maxDt);
            steps = std::max(1, std::min(steps, 100));
            float subDt = dt / steps;

            for (int s = 0; s < steps; s++) {
                attractors[i].step(subDt);
            }

            // Get raw normalized outputs (-1 to +1)
            float rawX = attractors[i].getNormX() / 5.0f;  // getNormX returns ±5V, convert to ±1
            float rawY = attractors[i].getNormY() / 5.0f;
            float rawZ = attractors[i].getNormZ() / 5.0f;

            // Apply smoothing
            smoothedX[i] += smoothCoeff * (rawX - smoothedX[i]);
            smoothedY[i] += smoothCoeff * (rawY - smoothedY[i]);
            smoothedZ[i] += smoothCoeff * (rawZ - smoothedZ[i]);

            // Scale to voltage
            bankOutputs[i][0] = scaleVoltage(smoothedX[i], voltageMode);
            bankOutputs[i][1] = scaleVoltage(smoothedY[i], voltageMode);
            bankOutputs[i][2] = scaleVoltage(smoothedZ[i], voltageMode);
            bankOutputs[i][3] = bankOutputs[i][0] + bankOutputs[i][1] + bankOutputs[i][2];
        }

        // Bank A outputs
        outputs[A_X_OUTPUT].setVoltage(bankOutputs[0][0]);
        outputs[A_Y_OUTPUT].setVoltage(bankOutputs[0][1]);
        outputs[A_Z_OUTPUT].setVoltage(bankOutputs[0][2]);
        outputs[A_SUM_OUTPUT].setVoltage(bankOutputs[0][3]);

        // Bank B outputs
        outputs[B_X_OUTPUT].setVoltage(bankOutputs[1][0]);
        outputs[B_Y_OUTPUT].setVoltage(bankOutputs[1][1]);
        outputs[B_Z_OUTPUT].setVoltage(bankOutputs[1][2]);
        outputs[B_SUM_OUTPUT].setVoltage(bankOutputs[1][3]);

        // Bank C outputs
        outputs[C_X_OUTPUT].setVoltage(bankOutputs[2][0]);
        outputs[C_Y_OUTPUT].setVoltage(bankOutputs[2][1]);
        outputs[C_Z_OUTPUT].setVoltage(bankOutputs[2][2]);
        outputs[C_SUM_OUTPUT].setVoltage(bankOutputs[2][3]);

        // Bank D outputs
        outputs[D_X_OUTPUT].setVoltage(bankOutputs[3][0]);
        outputs[D_Y_OUTPUT].setVoltage(bankOutputs[3][1]);
        outputs[D_Z_OUTPUT].setVoltage(bankOutputs[3][2]);
        outputs[D_SUM_OUTPUT].setVoltage(bankOutputs[3][3]);

        // Combined outputs (using smoothed normalized values for consistency)
        float combSum = bankOutputs[0][3] + bankOutputs[1][3] + bankOutputs[2][3] + bankOutputs[3][3];
        float combRect = std::abs(bankOutputs[0][3]) + std::abs(bankOutputs[1][3]) + std::abs(bankOutputs[2][3]) + std::abs(bankOutputs[3][3]);
        float combInv = -combSum;
        float combDist = 5.f - std::abs(bankOutputs[0][3]) - std::abs(bankOutputs[1][3]) - std::abs(bankOutputs[2][3]) - std::abs(bankOutputs[3][3]);

        outputs[COMB_SUM_OUTPUT].setVoltage(combSum);
        outputs[COMB_RECT_OUTPUT].setVoltage(combRect);
        outputs[COMB_INV_OUTPUT].setVoltage(combInv);
        outputs[COMB_DIST_OUTPUT].setVoltage(combDist);

        // Update trail history (downsample for display)
        trailCounter++;
        if (trailCounter >= (int)(args.sampleRate / 60.f)) { // ~60 fps
            trailCounter = 0;
            trailIndex = (trailIndex + 1) % MAX_TRAIL_LENGTH;

            // Store smoothed normalized positions for display (-1 to 1 range, clamped)
            for (int i = 0; i < 4; i++) {
                trailX[i][trailIndex] = clamp(smoothedX[i], -1.f, 1.f);
                trailY[i][trailIndex] = clamp(smoothedY[i], -1.f, 1.f);
                trailZ[i][trailIndex] = clamp(smoothedZ[i], -1.f, 1.f);
            }

            // Combined: use sum and rectified sum as x,y,z
            float normSum = (smoothedX[0] + smoothedY[0] + smoothedZ[0] +
                            smoothedX[1] + smoothedY[1] + smoothedZ[1] +
                            smoothedX[2] + smoothedY[2] + smoothedZ[2] +
                            smoothedX[3] + smoothedY[3] + smoothedZ[3]) / 12.f;
            float normRect = (std::abs(smoothedX[0]) + std::abs(smoothedY[0]) + std::abs(smoothedZ[0]) +
                             std::abs(smoothedX[1]) + std::abs(smoothedY[1]) + std::abs(smoothedZ[1]) +
                             std::abs(smoothedX[2]) + std::abs(smoothedY[2]) + std::abs(smoothedZ[2]) +
                             std::abs(smoothedX[3]) + std::abs(smoothedY[3]) + std::abs(smoothedZ[3])) / 12.f;
            combTrailX[trailIndex] = clamp(normSum, -1.f, 1.f);
            combTrailY[trailIndex] = clamp(normRect * 2.f - 1.f, -1.f, 1.f);
            combTrailZ[trailIndex] = clamp((smoothedZ[0] + smoothedZ[1] + smoothedZ[2] + smoothedZ[3]) / 4.f, -1.f, 1.f);
        }
    }

    json_t* dataToJson() override {
        json_t* rootJ = json_object();
        json_object_set_new(rootJ, "displayMode", json_integer(displayMode));
        json_object_set_new(rootJ, "display3D", json_boolean(display3D));
        json_object_set_new(rootJ, "displayStyle", json_integer(displayStyle));
        json_object_set_new(rootJ, "ajmanEnabled", json_boolean(ajmanEnabled));
        return rootJ;
    }

    void dataFromJson(json_t* rootJ) override {
        json_t* displayModeJ = json_object_get(rootJ, "displayMode");
        if (displayModeJ) {
            displayMode = json_integer_value(displayModeJ);
        }
        json_t* display3DJ = json_object_get(rootJ, "display3D");
        if (display3DJ) {
            display3D = json_boolean_value(display3DJ);
        }
        json_t* displayStyleJ = json_object_get(rootJ, "displayStyle");
        if (displayStyleJ) {
            displayStyle = json_integer_value(displayStyleJ);
        }
        json_t* ajmanEnabledJ = json_object_get(rootJ, "ajmanEnabled");
        if (ajmanEnabledJ) {
            ajmanEnabled = json_boolean_value(ajmanEnabledJ);
        }
    }
};


// Custom display widget for attractor visualization
struct AttractorDisplay : Widget {
    StrangeWeather* module = nullptr;
    float rotationTime = 0.f;
    int ajmanImage = -1;  // NanoVG image handle for Ajman

    // 3D projection: rotate point and return screen coordinates
    void project3D(float x, float y, float z, float angleX, float angleY, float& screenX, float& screenY, float& depth) {
        // Rotate around Y axis
        float rx = x * std::cos(angleY) - z * std::sin(angleY);
        float rz = x * std::sin(angleY) + z * std::cos(angleY);

        // Rotate around X axis
        float ry = y * std::cos(angleX) - rz * std::sin(angleX);
        float finalZ = y * std::sin(angleX) + rz * std::cos(angleX);

        // Simple perspective
        float perspective = 1.0f / (1.0f + finalZ * 0.15f);
        screenX = rx * perspective;
        screenY = ry * perspective;
        depth = finalZ;
    }

    void draw(const DrawArgs& args) override {
        // Update rotation time for 3D mode
        rotationTime += 1.f / 60.f;  // Approximate 60fps

        // Background
        nvgBeginPath(args.vg);
        nvgRect(args.vg, 0, 0, box.size.x, box.size.y);
        nvgFillColor(args.vg, nvgRGB(0x11, 0x11, 0x11));
        nvgFill(args.vg);

        // Clip drawing to display bounds
        nvgSave(args.vg);
        nvgScissor(args.vg, 0, 0, box.size.x, box.size.y);

        if (!module) {
            drawPreview(args);
            nvgRestore(args.vg);
            return;
        }

        int mode = module->displayMode;

        // Easter egg: Ajman mode (mode 6)
        if (mode == 6) {
            // Load image if not loaded yet
            if (ajmanImage == -1) {
                std::string path = asset::plugin(pluginInstance, "res/ajman.jpg");
                ajmanImage = nvgCreateImage(args.vg, path.c_str(), 0);
            }
            if (ajmanImage != -1) {
                // Draw image scaled to fit display with border
                float border = 4.f;
                NVGpaint imgPaint = nvgImagePattern(args.vg, border, border,
                    box.size.x - border * 2, box.size.y - border * 2, 0, ajmanImage, 1.0f);
                nvgBeginPath(args.vg);
                nvgRect(args.vg, border, border, box.size.x - border * 2, box.size.y - border * 2);
                nvgFillPaint(args.vg, imgPaint);
                nvgFill(args.vg);
            }
            nvgRestore(args.vg);
            return;
        }
        bool is3D = module->display3D;

        // Setup label style
        nvgFontSize(args.vg, 8);
        nvgFontFaceId(args.vg, APP->window->uiFont->handle);
        nvgFillColor(args.vg, nvgRGBA(0xff, 0xff, 0xff, 0x88));

        if (mode == 5) {
            // All four banks in quadrants (A, B, C, D)
            float w = box.size.x / 2.f;
            float h = box.size.y / 2.f;

            drawAttractor(args, 0, 0, 0, w, h, nvgRGB(0x00, 0xff, 0xaa), is3D);
            drawAttractor(args, 1, w, 0, w, h, nvgRGB(0xff, 0xaa, 0x00), is3D);
            drawAttractor(args, 2, 0, h, w, h, nvgRGB(0xaa, 0x00, 0xff), is3D);
            drawAttractor(args, 3, w, h, w, h, nvgRGB(0xff, 0x00, 0x66), is3D);

            // Grid lines
            nvgBeginPath(args.vg);
            nvgMoveTo(args.vg, w, 0);
            nvgLineTo(args.vg, w, box.size.y);
            nvgMoveTo(args.vg, 0, h);
            nvgLineTo(args.vg, box.size.x, h);
            nvgStrokeColor(args.vg, nvgRGBA(0x33, 0x33, 0x33, 0xff));
            nvgStrokeWidth(args.vg, 1.f);
            nvgStroke(args.vg);

            // Quadrant labels (bottom-left of each quadrant)
            nvgFillColor(args.vg, nvgRGBA(0xff, 0xff, 0xff, 0x88));
            nvgTextAlign(args.vg, NVG_ALIGN_LEFT | NVG_ALIGN_BOTTOM);
            nvgText(args.vg, 3, h - 3, "A", NULL);
            nvgText(args.vg, w + 3, h - 3, "B", NULL);
            nvgText(args.vg, 3, box.size.y - 3, "C", NULL);
            nvgText(args.vg, w + 3, box.size.y - 3, "D", NULL);
        }
        else if (mode == 4) {
            // Combined view
            drawCombined(args, 0, 0, box.size.x, box.size.y, nvgRGB(0xff, 0xff, 0xff), is3D);
            // Label
            nvgFillColor(args.vg, nvgRGBA(0xff, 0xff, 0xff, 0x88));
            nvgTextAlign(args.vg, NVG_ALIGN_RIGHT | NVG_ALIGN_TOP);
            nvgText(args.vg, box.size.x - 3, 3, "MIX", NULL);
        }
        else {
            // Single bank view (0=A, 1=B, 2=C, 3=D)
            NVGcolor colors[4] = {
                nvgRGB(0x00, 0xff, 0xaa),
                nvgRGB(0xff, 0xaa, 0x00),
                nvgRGB(0xaa, 0x00, 0xff),
                nvgRGB(0xff, 0x00, 0x66)
            };
            const char* labels[4] = {"A", "B", "C", "D"};
            drawAttractor(args, mode, 0, 0, box.size.x, box.size.y, colors[mode], is3D);
            // Label
            nvgFillColor(args.vg, nvgRGBA(0xff, 0xff, 0xff, 0x88));
            nvgTextAlign(args.vg, NVG_ALIGN_RIGHT | NVG_ALIGN_TOP);
            nvgText(args.vg, box.size.x - 3, 3, labels[mode], NULL);
        }

        // Show 3D indicator
        if (is3D) {
            nvgFillColor(args.vg, nvgRGBA(0xff, 0xff, 0xff, 0x88));
            nvgTextAlign(args.vg, NVG_ALIGN_LEFT | NVG_ALIGN_TOP);
            nvgText(args.vg, 3, 3, "3D", NULL);
        }

        // Restore clipping
        nvgRestore(args.vg);
    }

    void drawPreview(const DrawArgs& args) {
        nvgBeginPath(args.vg);
        float cx = box.size.x / 2.f;
        float cy = box.size.y / 2.f;
        float r = std::min(cx, cy) * 0.6f;

        for (int i = 0; i < 100; i++) {
            float t = i / 100.f * 2.f * M_PI;
            float x = cx + r * std::sin(t) * std::cos(t * 0.5f);
            float y = cy + r * std::cos(t) * 0.7f;
            if (i == 0)
                nvgMoveTo(args.vg, x, y);
            else
                nvgLineTo(args.vg, x, y);
        }
        nvgStrokeColor(args.vg, nvgRGBA(0x00, 0xff, 0xaa, 0x88));
        nvgStrokeWidth(args.vg, 1.f);
        nvgStroke(args.vg);
    }

    void drawAttractor(const DrawArgs& args, int bank, float ox, float oy, float w, float h, NVGcolor color, bool is3D) {
        if (!module) return;

        float cx = ox + w / 2.f;
        float cy = oy + h / 2.f;
        float scale = std::min(w, h) / 2.f * 0.75f;  // 75% to stay within bounds

        int idx = module->trailIndex;
        int style = module->displayStyle;

        // Rotation angles for 3D mode
        float angleX = rotationTime * 0.2f;
        float angleY = rotationTime * 0.33f;

        int trailLen = module->getTrailLength();

        if (style == 2) {
            // Scope mode: time-based waveforms (X, Y, Z stacked)
            float rowH = h / 3.f;
            float* data[3] = {module->trailX[bank], module->trailY[bank], module->trailZ[bank]};

            for (int row = 0; row < 3; row++) {
                float baseY = oy + rowH * row + rowH / 2.f;
                float waveScale = rowH * 0.4f;

                nvgBeginPath(args.vg);
                for (int i = 0; i < trailLen; i++) {
                    int i0 = (idx - trailLen + 1 + i + StrangeWeather::MAX_TRAIL_LENGTH) % StrangeWeather::MAX_TRAIL_LENGTH;
                    float xPos = ox + (float)i / trailLen * w;
                    float yPos = baseY + data[row][i0] * waveScale;
                    if (i == 0)
                        nvgMoveTo(args.vg, xPos, yPos);
                    else
                        nvgLineTo(args.vg, xPos, yPos);
                }
                nvgStrokeColor(args.vg, nvgRGBAf(color.r, color.g, color.b, 0.8f));
                nvgStrokeWidth(args.vg, 1.f);
                nvgStroke(args.vg);
            }
        }
        else if (style == 1 || is3D) {
            // Lissajous mode (style 1) or 3D mode: tiny dots
            for (int i = 0; i < trailLen; i++) {
                int i0 = (idx - i + StrangeWeather::MAX_TRAIL_LENGTH) % StrangeWeather::MAX_TRAIL_LENGTH;

                float x, y;
                if (is3D) {
                    float sx, sy, depth;
                    project3D(module->trailX[bank][i0], module->trailY[bank][i0], module->trailZ[bank][i0],
                             angleX, angleY, sx, sy, depth);
                    x = cx + sx * scale;
                    y = cy + sy * scale;
                } else {
                    x = cx + module->trailX[bank][i0] * scale;
                    y = cy + module->trailY[bank][i0] * scale;
                }

                // Age-based fade
                float age = (float)i / trailLen;
                float brightness = 1.f - age;
                brightness = brightness * brightness;

                // Tiny dots for Lissajous phosphor look
                float dotSize = (style == 1) ? 0.8f : 1.5f;

                nvgBeginPath(args.vg);
                nvgCircle(args.vg, x, y, dotSize);
                nvgFillColor(args.vg, nvgRGBAf(color.r, color.g, color.b, brightness * 0.9f));
                nvgFill(args.vg);
            }
        }
        else {
            // Trace mode (style 0): lines
            for (int i = 0; i < trailLen - 1; i++) {
                int i0 = (idx - i + StrangeWeather::MAX_TRAIL_LENGTH) % StrangeWeather::MAX_TRAIL_LENGTH;
                int i1 = (idx - i - 1 + StrangeWeather::MAX_TRAIL_LENGTH) % StrangeWeather::MAX_TRAIL_LENGTH;

                float alpha = 1.f - (float)i / trailLen;
                alpha = alpha * alpha * 0.8f;

                float x0 = cx + module->trailX[bank][i0] * scale;
                float y0 = cy + module->trailY[bank][i0] * scale;
                float x1 = cx + module->trailX[bank][i1] * scale;
                float y1 = cy + module->trailY[bank][i1] * scale;

                nvgBeginPath(args.vg);
                nvgMoveTo(args.vg, x0, y0);
                nvgLineTo(args.vg, x1, y1);
                nvgStrokeColor(args.vg, nvgRGBAf(color.r, color.g, color.b, alpha));
                nvgStrokeWidth(args.vg, 1.f + alpha);
                nvgStroke(args.vg);
            }

            // Current position dot
            float x = cx + module->trailX[bank][idx] * scale;
            float y = cy + module->trailY[bank][idx] * scale;
            nvgBeginPath(args.vg);
            nvgCircle(args.vg, x, y, 2.f);
            nvgFillColor(args.vg, color);
            nvgFill(args.vg);
        }
    }
    
    void drawCombined(const DrawArgs& args, float ox, float oy, float w, float h, NVGcolor color, bool is3D) {
        if (!module) return;

        float cx = ox + w / 2.f;
        float cy = oy + h / 2.f;
        float scale = std::min(w, h) / 2.f * 0.75f;  // 75% to stay within bounds

        int idx = module->trailIndex;
        int style = module->displayStyle;

        // Rotation angles for 3D mode (slower, more contemplative)
        float angleX = rotationTime * 0.2f;
        float angleY = rotationTime * 0.33f;

        int trailLen = module->getTrailLength();

        if (style == 2) {
            // Scope mode: time-based waveforms (X, Y, Z stacked)
            float rowH = h / 3.f;
            float* data[3] = {module->combTrailX, module->combTrailY, module->combTrailZ};

            for (int row = 0; row < 3; row++) {
                float baseY = oy + rowH * row + rowH / 2.f;
                float waveScale = rowH * 0.4f;

                nvgBeginPath(args.vg);
                for (int i = 0; i < trailLen; i++) {
                    int i0 = (idx - trailLen + 1 + i + StrangeWeather::MAX_TRAIL_LENGTH) % StrangeWeather::MAX_TRAIL_LENGTH;
                    float xPos = ox + (float)i / trailLen * w;
                    float yPos = baseY + data[row][i0] * waveScale;
                    if (i == 0)
                        nvgMoveTo(args.vg, xPos, yPos);
                    else
                        nvgLineTo(args.vg, xPos, yPos);
                }
                nvgStrokeColor(args.vg, nvgRGBAf(color.r, color.g, color.b, 0.8f));
                nvgStrokeWidth(args.vg, 1.f);
                nvgStroke(args.vg);
            }
        }
        else if (style == 1 || is3D) {
            // Lissajous mode (style 1) or 3D mode: tiny dots
            for (int i = 0; i < trailLen; i++) {
                int i0 = (idx - i + StrangeWeather::MAX_TRAIL_LENGTH) % StrangeWeather::MAX_TRAIL_LENGTH;

                float x, y;
                if (is3D) {
                    float sx, sy, depth;
                    project3D(module->combTrailX[i0], module->combTrailY[i0], module->combTrailZ[i0],
                             angleX, angleY, sx, sy, depth);
                    x = cx + sx * scale;
                    y = cy + sy * scale;
                } else {
                    x = cx + module->combTrailX[i0] * scale;
                    y = cy + module->combTrailY[i0] * scale;
                }

                // Age-based fade
                float age = (float)i / trailLen;
                float brightness = 1.f - age;
                brightness = brightness * brightness;

                // Tiny dots for Lissajous phosphor look
                float dotSize = (style == 1) ? 0.8f : 1.5f;

                nvgBeginPath(args.vg);
                nvgCircle(args.vg, x, y, dotSize);
                nvgFillColor(args.vg, nvgRGBAf(color.r, color.g, color.b, brightness * 0.9f));
                nvgFill(args.vg);
            }
        }
        else {
            // Trace mode (style 0): lines
            for (int i = 0; i < trailLen - 1; i++) {
                int i0 = (idx - i + StrangeWeather::MAX_TRAIL_LENGTH) % StrangeWeather::MAX_TRAIL_LENGTH;
                int i1 = (idx - i - 1 + StrangeWeather::MAX_TRAIL_LENGTH) % StrangeWeather::MAX_TRAIL_LENGTH;

                float alpha = 1.f - (float)i / trailLen;
                alpha = alpha * alpha * 0.8f;

                float x0 = cx + module->combTrailX[i0] * scale;
                float y0 = cy + module->combTrailY[i0] * scale;
                float x1 = cx + module->combTrailX[i1] * scale;
                float y1 = cy + module->combTrailY[i1] * scale;

                nvgBeginPath(args.vg);
                nvgMoveTo(args.vg, x0, y0);
                nvgLineTo(args.vg, x1, y1);
                nvgStrokeColor(args.vg, nvgRGBAf(color.r, color.g, color.b, alpha));
                nvgStrokeWidth(args.vg, 1.f + alpha);
                nvgStroke(args.vg);
            }

            // Current position dot
            float x = cx + module->combTrailX[idx] * scale;
            float y = cy + module->combTrailY[idx] * scale;
            nvgBeginPath(args.vg);
            nvgCircle(args.vg, x, y, 2.f);
            nvgFillColor(args.vg, color);
            nvgFill(args.vg);
        }
    }
};


// Cycle button - triggers display mode change on click
struct CycleButton : app::SvgSwitch {
    StrangeWeather* swModule = nullptr;

    CycleButton() {
        momentary = true;
        addFrame(Svg::load(asset::system("res/ComponentLibrary/TL1105_0.svg")));
        addFrame(Svg::load(asset::system("res/ComponentLibrary/TL1105_1.svg")));
    }

    void onDragStart(const DragStartEvent& e) override {
        SvgSwitch::onDragStart(e);
        if (swModule) {
            swModule->cycleDisplay();
        }
    }
};

// Mode button - cycles through display styles (Trace/Lissajous/Scope)
struct ModeButton : app::SvgSwitch {
    StrangeWeather* swModule = nullptr;

    ModeButton() {
        momentary = true;
        addFrame(Svg::load(asset::system("res/ComponentLibrary/TL1105_0.svg")));
        addFrame(Svg::load(asset::system("res/ComponentLibrary/TL1105_1.svg")));
    }

    void onDragStart(const DragStartEvent& e) override {
        SvgSwitch::onDragStart(e);
        if (swModule) {
            swModule->cycleDisplayStyle();
        }
    }
};

// 3D toggle button - toggles 3D display mode
struct Toggle3DButton : app::SvgSwitch {
    StrangeWeather* swModule = nullptr;

    Toggle3DButton() {
        momentary = true;
        addFrame(Svg::load(asset::system("res/ComponentLibrary/TL1105_0.svg")));
        addFrame(Svg::load(asset::system("res/ComponentLibrary/TL1105_1.svg")));
    }

    void onDragStart(const DragStartEvent& e) override {
        SvgSwitch::onDragStart(e);
        if (swModule) {
            swModule->display3D = !swModule->display3D;
        }
    }
};

// 4-position vertical toggle switch
struct CKSSFour : app::SvgSwitch {
    CKSSFour() {
        shadow->opacity = 0.0;
        addFrame(Svg::load(asset::plugin(pluginInstance, "res/CKSSFour_0.svg")));
        addFrame(Svg::load(asset::plugin(pluginInstance, "res/CKSSFour_1.svg")));
        addFrame(Svg::load(asset::plugin(pluginInstance, "res/CKSSFour_2.svg")));
        addFrame(Svg::load(asset::plugin(pluginInstance, "res/CKSSFour_3.svg")));
    }
};

// Davies-style knob for rate control (like Synthesis Technology)
struct DaviesKnob : Davies1900hBlackKnob {
    DaviesKnob() {}
};

// Small Davies knob for chaos
struct SmallDaviesKnob : Davies1900hBlackKnob {
    SmallDaviesKnob() {
        // Scale down
        box.size = mm2px(Vec(8.0, 8.0));
    }
};

// Custom widget to draw panel labels using NanoVG
struct PanelLabels : Widget {
    void draw(const DrawArgs& args) override {
        nvgFontSize(args.vg, 10);
        nvgFontFaceId(args.vg, APP->window->uiFont->handle);
        nvgFillColor(args.vg, nvgRGB(0x22, 0x22, 0x22));

        // Module title (centered on 40HP panel = 101.6mm center)
        nvgFontSize(args.vg, 12);
        nvgTextAlign(args.vg, NVG_ALIGN_CENTER | NVG_ALIGN_TOP);
        nvgText(args.vg, mm2px(101.6), mm2px(2), "STRANGE WEATHER", NULL);

        // Column headers
        nvgFontSize(args.vg, 8);
        nvgFillColor(args.vg, nvgRGB(0x44, 0x44, 0x44));
        float headerY = mm2px(10);
        nvgText(args.vg, mm2px(85), headerY, "RATE", NULL);
        nvgText(args.vg, mm2px(100), headerY, "RNG", NULL);
        nvgText(args.vg, mm2px(112), headerY, "SHAPE", NULL);
        nvgText(args.vg, mm2px(124), headerY, "VOLT", NULL);
        nvgText(args.vg, mm2px(139), headerY, "CHAOS", NULL);
        nvgText(args.vg, mm2px(158), headerY, "x", NULL);
        nvgText(args.vg, mm2px(170), headerY, "y", NULL);
        nvgText(args.vg, mm2px(182), headerY, "z", NULL);
        nvgText(args.vg, mm2px(194), headerY, "SUM", NULL);

        // Bank labels - single letters, to left of controls
        nvgFontSize(args.vg, 10);
        nvgTextAlign(args.vg, NVG_ALIGN_RIGHT | NVG_ALIGN_MIDDLE);
        nvgFillColor(args.vg, nvgRGB(0x00, 0x99, 0x66));
        nvgText(args.vg, mm2px(77), mm2px(24), "A", NULL);
        nvgFillColor(args.vg, nvgRGB(0xcc, 0x88, 0x00));
        nvgText(args.vg, mm2px(77), mm2px(48), "B", NULL);
        nvgFillColor(args.vg, nvgRGB(0x88, 0x00, 0xcc));
        nvgText(args.vg, mm2px(77), mm2px(72), "C", NULL);
        nvgFillColor(args.vg, nvgRGB(0xcc, 0x00, 0x44));
        nvgText(args.vg, mm2px(77), mm2px(96), "D", NULL);

        // CYCLE, MODE, 3D, and TRAIL labels (below display, centered on 35mm)
        nvgFontSize(args.vg, 8);
        nvgTextAlign(args.vg, NVG_ALIGN_CENTER | NVG_ALIGN_TOP);
        nvgFillColor(args.vg, nvgRGB(0x44, 0x44, 0x44));
        nvgText(args.vg, mm2px(17.0), mm2px(82), "CYCLE", NULL);
        nvgText(args.vg, mm2px(29.0), mm2px(82), "MODE", NULL);
        nvgText(args.vg, mm2px(41.0), mm2px(82), "3D", NULL);
        nvgText(args.vg, mm2px(53.0), mm2px(82), "TRAIL", NULL);

        // Attractor type descriptions (under display, in line with COMBINED)
        nvgFontSize(args.vg, 7);
        nvgTextAlign(args.vg, NVG_ALIGN_LEFT | NVG_ALIGN_TOP);
        nvgFillColor(args.vg, nvgRGB(0x55, 0x55, 0x55));
        nvgText(args.vg, mm2px(5), mm2px(90), "LORENZ: Smooth two-lobed butterfly", NULL);
        nvgText(args.vg, mm2px(5), mm2px(95), "ROSSLER: Asymmetric spiral, large excursions", NULL);
        nvgText(args.vg, mm2px(5), mm2px(100), "THOMAS: Cyclically symmetric, smooth rolling", NULL);
        nvgText(args.vg, mm2px(5), mm2px(105), "HALVORSEN: Sculptural, sharp transitions", NULL);

        // Combined section - label under RATE column, inline with jacks
        nvgFontSize(args.vg, 10);
        nvgTextAlign(args.vg, NVG_ALIGN_CENTER | NVG_ALIGN_MIDDLE);
        nvgFillColor(args.vg, nvgRGB(0x99, 0x66, 0x00));
        nvgText(args.vg, mm2px(85), mm2px(120), "COMBINED", NULL);

        // Combined output labels
        nvgFontSize(args.vg, 8);
        nvgTextAlign(args.vg, NVG_ALIGN_CENTER | NVG_ALIGN_TOP);
        nvgFillColor(args.vg, nvgRGB(0x44, 0x44, 0x44));
        float combY = mm2px(113);
        nvgText(args.vg, mm2px(115), combY, "SUM", NULL);
        nvgText(args.vg, mm2px(135), combY, "RECT", NULL);
        nvgText(args.vg, mm2px(155), combY, "INV", NULL);
        nvgText(args.vg, mm2px(175), combY, "DIST", NULL);
    }
};


struct StrangeWeatherWidget : ModuleWidget {
    AttractorDisplay* display = nullptr;
    
    StrangeWeatherWidget(StrangeWeather* module) {
        setModule(module);
        setPanel(createPanel(asset::plugin(pluginInstance, "res/StrangeWeather.svg")));
        
        // Screws
        addChild(createWidget<ScrewSilver>(Vec(RACK_GRID_WIDTH, 0)));
        addChild(createWidget<ScrewSilver>(Vec(box.size.x - 2 * RACK_GRID_WIDTH, 0)));
        addChild(createWidget<ScrewSilver>(Vec(RACK_GRID_WIDTH, RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));
        addChild(createWidget<ScrewSilver>(Vec(box.size.x - 2 * RACK_GRID_WIDTH, RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));

        // Panel labels (drawn with NanoVG)
        PanelLabels* labels = new PanelLabels();
        labels->box.pos = Vec(0, 0);
        labels->box.size = box.size;
        addChild(labels);

        // Display (large 60x60mm square)
        display = new AttractorDisplay();
        display->box.pos = mm2px(Vec(5.0, 12.0));
        display->box.size = mm2px(Vec(60.0, 60.0));
        display->module = module;
        addChild(display);

        // Cycle button (below display, centered on 35mm with 12mm spacing)
        CycleButton* cycleBtn = new CycleButton();
        cycleBtn->box.pos = mm2px(Vec(17.0, 77.0)).minus(cycleBtn->box.size.div(2));
        cycleBtn->swModule = module;
        addChild(cycleBtn);

        // Mode button (cycles display style)
        ModeButton* modeBtn = new ModeButton();
        modeBtn->box.pos = mm2px(Vec(29.0, 77.0)).minus(modeBtn->box.size.div(2));
        modeBtn->swModule = module;
        addChild(modeBtn);

        // 3D toggle button
        Toggle3DButton* toggle3DBtn = new Toggle3DButton();
        toggle3DBtn->box.pos = mm2px(Vec(41.0, 77.0)).minus(toggle3DBtn->box.size.div(2));
        toggle3DBtn->swModule = module;
        addChild(toggle3DBtn);

        // Trail length knob
        addParam(createParamCentered<Trimpot>(mm2px(Vec(53.0, 77.0)), module, StrangeWeather::TRAIL_PARAM));

        // Layout: Each bank has Rate knob, Range toggle (3-pos), Shape toggle (4-pos),
        //         Voltage toggle (4-pos), Chaos knob, then 4 outputs
        // All x positions shifted +30mm for 36HP panel

        // Bank A controls and outputs (y = 24mm center)
        // 40HP panel = 203.2mm wide. Controls start at 85mm to leave room for labels
        float yA = 24.0;
        addParam(createParamCentered<DaviesKnob>(mm2px(Vec(85.0, yA)), module, StrangeWeather::RATE_A_PARAM));
        addParam(createParamCentered<CKSSThree>(mm2px(Vec(100.0, yA)), module, StrangeWeather::RANGE_A_PARAM));
        addParam(createParamCentered<CKSSFour>(mm2px(Vec(112.0, yA)), module, StrangeWeather::SHAPE_A_PARAM));
        addParam(createParamCentered<CKSSFour>(mm2px(Vec(124.0, yA)), module, StrangeWeather::VOLTAGE_A_PARAM));
        addParam(createParamCentered<DaviesKnob>(mm2px(Vec(139.0, yA)), module, StrangeWeather::CHAOS_A_PARAM));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(158.0, yA)), module, StrangeWeather::A_X_OUTPUT));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(170.0, yA)), module, StrangeWeather::A_Y_OUTPUT));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(182.0, yA)), module, StrangeWeather::A_Z_OUTPUT));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(194.0, yA)), module, StrangeWeather::A_SUM_OUTPUT));

        // Bank B controls and outputs (y = 48mm center)
        float yB = 48.0;
        addParam(createParamCentered<DaviesKnob>(mm2px(Vec(85.0, yB)), module, StrangeWeather::RATE_B_PARAM));
        addParam(createParamCentered<CKSSThree>(mm2px(Vec(100.0, yB)), module, StrangeWeather::RANGE_B_PARAM));
        addParam(createParamCentered<CKSSFour>(mm2px(Vec(112.0, yB)), module, StrangeWeather::SHAPE_B_PARAM));
        addParam(createParamCentered<CKSSFour>(mm2px(Vec(124.0, yB)), module, StrangeWeather::VOLTAGE_B_PARAM));
        addParam(createParamCentered<DaviesKnob>(mm2px(Vec(139.0, yB)), module, StrangeWeather::CHAOS_B_PARAM));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(158.0, yB)), module, StrangeWeather::B_X_OUTPUT));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(170.0, yB)), module, StrangeWeather::B_Y_OUTPUT));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(182.0, yB)), module, StrangeWeather::B_Z_OUTPUT));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(194.0, yB)), module, StrangeWeather::B_SUM_OUTPUT));

        // Bank C controls and outputs (y = 72mm center)
        float yC = 72.0;
        addParam(createParamCentered<DaviesKnob>(mm2px(Vec(85.0, yC)), module, StrangeWeather::RATE_C_PARAM));
        addParam(createParamCentered<CKSSThree>(mm2px(Vec(100.0, yC)), module, StrangeWeather::RANGE_C_PARAM));
        addParam(createParamCentered<CKSSFour>(mm2px(Vec(112.0, yC)), module, StrangeWeather::SHAPE_C_PARAM));
        addParam(createParamCentered<CKSSFour>(mm2px(Vec(124.0, yC)), module, StrangeWeather::VOLTAGE_C_PARAM));
        addParam(createParamCentered<DaviesKnob>(mm2px(Vec(139.0, yC)), module, StrangeWeather::CHAOS_C_PARAM));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(158.0, yC)), module, StrangeWeather::C_X_OUTPUT));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(170.0, yC)), module, StrangeWeather::C_Y_OUTPUT));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(182.0, yC)), module, StrangeWeather::C_Z_OUTPUT));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(194.0, yC)), module, StrangeWeather::C_SUM_OUTPUT));

        // Bank D controls and outputs (y = 96mm center)
        float yD = 96.0;
        addParam(createParamCentered<DaviesKnob>(mm2px(Vec(85.0, yD)), module, StrangeWeather::RATE_D_PARAM));
        addParam(createParamCentered<CKSSThree>(mm2px(Vec(100.0, yD)), module, StrangeWeather::RANGE_D_PARAM));
        addParam(createParamCentered<CKSSFour>(mm2px(Vec(112.0, yD)), module, StrangeWeather::SHAPE_D_PARAM));
        addParam(createParamCentered<CKSSFour>(mm2px(Vec(124.0, yD)), module, StrangeWeather::VOLTAGE_D_PARAM));
        addParam(createParamCentered<DaviesKnob>(mm2px(Vec(139.0, yD)), module, StrangeWeather::CHAOS_D_PARAM));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(158.0, yD)), module, StrangeWeather::D_X_OUTPUT));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(170.0, yD)), module, StrangeWeather::D_Y_OUTPUT));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(182.0, yD)), module, StrangeWeather::D_Z_OUTPUT));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(194.0, yD)), module, StrangeWeather::D_SUM_OUTPUT));

        // Combined outputs (y = 120mm center, shifted right)
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(115.0, 120.0)), module, StrangeWeather::COMB_SUM_OUTPUT));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(135.0, 120.0)), module, StrangeWeather::COMB_RECT_OUTPUT));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(155.0, 120.0)), module, StrangeWeather::COMB_INV_OUTPUT));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(175.0, 120.0)), module, StrangeWeather::COMB_DIST_OUTPUT));
    }
    
    void step() override {
        ModuleWidget::step();
        
        // Handle cycle button (hacky but works for now)
        // A better approach would be a custom button widget
    }
    
    void appendContextMenu(Menu* menu) override {
        StrangeWeather* module = dynamic_cast<StrangeWeather*>(this->module);
        if (!module) return;

        menu->addChild(new MenuSeparator());
        menu->addChild(createMenuLabel("Display"));

        menu->addChild(createIndexSubmenuItem("View",
            {"Bank A", "Bank B", "Bank C", "Bank D", "Combined", "All"},
            [=]() { return module->displayMode; },
            [=](int mode) { module->displayMode = mode; }
        ));

        menu->addChild(createBoolPtrMenuItem("3D Rotation", "", &module->display3D));

        menu->addChild(new MenuSeparator());
        menu->addChild(createBoolPtrMenuItem("Ajman", "", &module->ajmanEnabled));
    }
};


Model* modelStrangeWeather = createModel<StrangeWeather, StrangeWeatherWidget>("StrangeWeather");
