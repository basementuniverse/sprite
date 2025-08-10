import { vec2 } from '@basementuniverse/vec';
export type SpriteOptionsData = Partial<Omit<SpriteOptions, 'image' | 'preRender' | 'postRender' | 'debug'>> & {
    imageName?: string;
    animations?: {
        [name: string]: {
            [direction: string]: SpriteAnimationOptionsData;
        };
    };
};
export type SpriteAnimationOptionsData = Omit<SpriteAnimationOptions, 'images'> & {
    imageNames?: string[];
};
export type SpriteOptions = {
    /**
     * The position of the sprite
     *
     * Defaults to (0, 0)
     */
    position?: vec2;
    /**
     * The base size of the sprite
     *
     * If omitted, use the base image size, or fall back to the size of the
     * first image found in available directions for the default animation
     *
     * If a size still can't be found, default to (0, 0)
     */
    size?: vec2;
    /**
     * Origin offset from top-left corner, used for rotation and scaling
     *
     * If omitted, this will be placed in the center of the sprite (based on
     * size, noting that size will be calculated first)
     *
     * @see SpriteOptions.size
     */
    origin?: vec2;
    /**
     * The scale factor of the sprite
     *
     * Default is 1
     */
    scale?: number;
    /**
     * The sprite rotation, measured in radians
     *
     * Default is 0
     */
    rotation?: number;
    /**
     * An array of valid direction names
     *
     * By default a sprite will have one available direction: 'default'
     */
    directions: string[];
    /**
     * The initial direction of the sprite
     *
     * Default is 'default'
     */
    defaultDirection: string;
    /**
     * An optional base image
     *
     * If an animation frame doesn't exist (for example an animation has been
     * configured with 5 frames but only 3 images exist in the animation's
     * images array), we fall back to this image
     *
     * If an animation frame image can't be found and this base image doesn't
     * exist, we can still render the sprite - it will just have an empty image
     *
     * (this can be useful for sprites that act purely as hosts for attachment
     * points, for example)
     */
    image?: HTMLImageElement | HTMLCanvasElement;
    /**
     * A dictionary of animations
     *
     * Each animation can have multiple variants for different directions
     *
     * Use '*' as a direction name to indicate that a variant is available for
     * all directions and can be used as a fallback if the current direction
     * can't be found
     */
    animations: {
        [name: string]: {
            [direction: string]: SpriteAnimationOptions;
        };
    };
    /**
     * The initial animation for the sprite
     */
    defaultAnimation: string;
    /**
     * A list of attachment points
     */
    attachmentPoints?: SpriteAttachmentPointOptions[];
    /**
     * Optional hook called before rendering the sprite image
     */
    preRender?: (context: CanvasRenderingContext2D, sprite: Sprite) => void;
    /**
     * Optional hook called after rendering the sprite image
     */
    postRender?: (context: CanvasRenderingContext2D, sprite: Sprite) => void;
    /**
     * Optional debug options
     *
     * Can be a boolean value (in which case all sub-options will be set to the
     * same value), or an object allowing specific debug options to be enabled
     * individually
     */
    debug?: Partial<SpriteDebugOptions> | boolean;
};
type SpriteDebugOptions = {
    showSpriteTransforms: boolean;
    showSpriteBoundingBox: boolean;
    showAttachmentPoints: boolean;
};
export declare enum SpriteAnimationRepeatMode {
    /**
     * Loop this animation indefinitely
     */
    Repeat = "repeat",
    /**
     * Play once and then stop on the last frame
     */
    PlayOnceAndStop = "play-once-and-stop",
    /**
     * Play once and then reset back to the first frame
     */
    PlayOnceAndReset = "play-once-and-reset"
}
export type SpriteAnimationOptions = {
    /**
     * The name of this animation
     */
    name: string;
    /**
     * The number of frames in this animation
     *
     * Default is 1
     */
    frameCount?: number;
    /**
     * The number of frames per second when playing this animation
     *
     * Default is 1
     */
    frameRate?: number;
    /**
     * The repeat mode for this animation
     *
     * Default is SpriteAnimationRepeatMode.Repeat
     */
    mode?: SpriteAnimationRepeatMode;
    /**
     * A list of images to use for each frame
     *
     * If this list is longer than frameCount, some images won't be used
     *
     * If this list if shorter than frameCount, we fall back to the base image
     * or don't show an image for some frames
     */
    images?: (HTMLImageElement | HTMLCanvasElement)[];
    /**
     * Keyframes for attachment points
     */
    attachmentPointKeyframes?: {
        [attachmentPointName: string]: SpriteAttachmentPointKeyframe[];
    };
};
export type SpriteAttachmentPointOptions = {
    /**
     * The name of this attachment point
     */
    name: string;
    /**
     * The position offset of this attachment point, measured from the origin
     * position of this sprite
     */
    offset: vec2;
};
export type SpriteAttachmentPointKeyframe = {
    /**
     * The frame index (0-based) for this keyframe
     *
     * These should be defined in ascending frame order, however there can
     * be gaps
     *
     * Values will be interpolated between keyframes
     */
    frame: number;
    /**
     * The attachment point offset for this keyframe
     */
    offset: vec2;
};
export type SpriteAttachmentPointMap = Record<string, vec2>;
export declare function isSpriteOptionsData(value: unknown): value is SpriteOptionsData;
export declare class Sprite {
    private static readonly DEFAULT_OPTIONS;
    private static readonly DEFAULT_ANIMATION_OPTIONS;
    private static readonly DEBUG_BOUNDING_BOX_COLOUR;
    private static readonly DEBUG_BOUNDING_BOX_LINE_WIDTH;
    private static readonly DEBUG_TRANSFORMS_COLOUR_X;
    private static readonly DEBUG_TRANSFORMS_COLOUR_Y;
    private static readonly DEBUG_TRANSFORMS_LINE_WIDTH;
    private static readonly DEBUG_TRANSFORMS_SIZE;
    private static readonly DEBUG_ATTACHMENT_POINT_COLOUR;
    private static readonly DEBUG_ATTACHMENT_POINT_LINE_WIDTH;
    private static readonly DEBUG_ATTACHMENT_POINT_SIZE;
    private options;
    position: vec2;
    size: vec2;
    origin: vec2;
    scale: number;
    rotation: number;
    private _direction;
    private _animation;
    private currentAnimationOptions;
    private currentAnimationState;
    private currentImage;
    private currentAttachmentPoints;
    constructor(options?: Partial<SpriteOptions>);
    get direction(): string;
    set direction(value: string);
    get animation(): string;
    set animation(value: string);
    playAnimation(): void;
    pauseAnimation(): void;
    resetAnimation(): void;
    getAttachmentPoint(name: string): vec2 | null;
    update(dt: number): void;
    private updateAnimationOptions;
    private updateAnimationState;
    private updateImage;
    private updateAttachmentPoints;
    private findPreviousKeyframe;
    draw(context: CanvasRenderingContext2D): void;
    private drawTransformsMarker;
    private drawCross;
}
/**
 * Content Manager Processor wrapper which converts SpriteOptionsData into
 * SpriteOptions
 *
 * @see https://www.npmjs.com/package/@basementuniverse/content-manager
 */
export declare function spriteOptionsContentProcessor(content: Record<string, {
    name: string;
    type: string;
    content: any;
    status: string;
}>, data: {
    name: string;
    type: string;
    content: any;
    status: string;
}): Promise<void>;
export {};
