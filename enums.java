/**
 * Enum types defined for lineardynamicalsystems.pde
 *
 * Since Processing does not fully support enum types, they must be defined
 * here.
 */


/**
 * Enum type for current mode of the app
 *    Menu
 *    Help
 *    Configuration
 *    Main (where the animation occurs)
 */
enum app_mode {
    MENU,
    HELP,
    CONFIG,
    MAIN
}

/**
 * Enum type for current field of the configuration option
 *
 * In the configuration mode, the user can change the follow properties.
 *    Dimension
 *    Matrix
 *    Vector
 *    Time Step Size
 *    Total Nodes
 *    Node Spread Ratio
 *    Continue
 * The Continue field acts as a jump to the animation; it is not an actual
 * editable field
 */
enum config_field {
    DIM,
    MAT,
    VEC,
    DT,
    NTOTAL,
    NSPREAD,
    CONT
}