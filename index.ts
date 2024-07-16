import { vec } from '@basementuniverse/vec';

export type SpriteOptionsData = Partial<Omit<
  SpriteOptions,
  | 'image'
  | 'preRender'
  | 'postRender'
  | 'debug'
>> & {
  imageName?: string;
  animations?: {
    [name: string]: {
      [direction: string]: SpriteAnimationOptionsData;
    };
  };
};

export type SpriteAnimationOptionsData = Omit<
  SpriteAnimationOptions, 'images'
> & {
  imageNames?: string[];
};

export type SpriteOptions = {
  /**
   * The position of the sprite
   *
   * Defaults to (0, 0)
   */
  position?: vec;

  /**
   * The base size of the sprite
   *
   * If omitted, use the base image size, or fall back to the size of the
   * first image found in available directions for the default animation
   *
   * If a size still can't be found, default to (0, 0)
   */
  size?: vec;

  /**
   * Origin offset from top-left corner, used for rotation and scaling
   *
   * If omitted, this will be placed in the center of the sprite (based on
   * size, noting that size will be calculated first)
   *
   * @see SpriteOptions.size
   */
  origin?: vec;

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
  preRender?: (
    context: CanvasRenderingContext2D,
    sprite: Sprite
  ) => void;

  /**
   * Optional hook called after rendering the sprite image
   */
  postRender?: (
    context: CanvasRenderingContext2D,
    sprite: Sprite
  ) => void;

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

export enum SpriteAnimationRepeatMode {
  /**
   * Loop this animation indefinitely
   */
  Repeat = 0,

  /**
   * Play once and then stop on the last frame
   */
  PlayOnceAndStop,

  /**
   * Play once and then reset back to the first frame
   */
  PlayOnceAndReset,
}

export type SpriteAnimationOptions = {
  /**
   * The name of this animation
   */
  name: string;

  /**
   * The number of frames in this animation
   */
  frameCount: number;

  /**
   * The number of frames per second when playing this animation
   */
  frameRate: number;

  /**
   * The repeat mode for this animation
   */
  mode: SpriteAnimationRepeatMode;

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

type SpriteAnimationState = {
  /**
   * This animation is currently playing
   */
  playing: boolean;

  /**
   * The current frame number that this animation is on
   */
  currentFrame: number;

  /**
   * The amount of time elapsed since the current frame began
   */
  currentFrameTime: number;
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
  offset: vec;
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
  offset: vec;
};

export type SpriteAttachmentPointMap = Record<string, vec>;

export class Sprite {
  private static readonly DEFAULT_OPTIONS: SpriteOptions = {
    directions: ['default'],
    defaultDirection: 'default',
    animations: {
      default: {
        '*': {
          name: 'default',
          frameCount: 1,
          frameRate: 1,
          mode: SpriteAnimationRepeatMode.PlayOnceAndStop,
        },
      },
    },
    defaultAnimation: 'default',
  };

  private static readonly DEBUG_BOUNDING_BOX_COLOUR = 'green';
  private static readonly DEBUG_BOUNDING_BOX_LINE_WIDTH = 2;

  private static readonly DEBUG_TRANSFORMS_COLOUR_X = 'red';
  private static readonly DEBUG_TRANSFORMS_COLOUR_Y = 'orange';
  private static readonly DEBUG_TRANSFORMS_LINE_WIDTH = 1;
  private static readonly DEBUG_TRANSFORMS_SIZE = 10;

  private static readonly DEBUG_ATTACHMENT_POINT_COLOUR = 'blue';
  private static readonly DEBUG_ATTACHMENT_POINT_LINE_WIDTH = 2;
  private static readonly DEBUG_ATTACHMENT_POINT_SIZE = 5;

  private options: SpriteOptions & {
    debug: Required<SpriteDebugOptions>;
  };

  public position: vec = vec();
  public size: vec = vec();

  public origin: vec = vec();
  public scale: number = 1;
  public rotation: number = 0;

  private _direction: string;
  private _animation: string;

  private currentAnimationOptions: SpriteAnimationOptions | null = null;
  private currentAnimationState: SpriteAnimationState | null = null;
  private currentImage: HTMLImageElement | HTMLCanvasElement | null = null;
  private currentAttachmentPoints: SpriteAttachmentPointMap | null = null;

  public constructor(options?: Partial<SpriteOptions>) {
    const actualOptions = Object.assign(
      {},
      Sprite.DEFAULT_OPTIONS,
      options ?? {}
    );

    if (!actualOptions.debug || actualOptions.debug === true) {
      actualOptions.debug = {
        showSpriteTransforms: !!actualOptions.debug,
        showSpriteBoundingBox: !!actualOptions.debug,
        showAttachmentPoints: !!actualOptions.debug,
      };
    }

    this.options = actualOptions as typeof this.options;

    if (this.options.position) {
      this.position = vec.cpy(this.options.position);
    }

    if (this.options.size) {
      this.size = vec.cpy(this.options.size);
    } else {
      // Default to the size of the base image if one exists
      if (this.options.image) {
        this.size = vec(
          this.options.image.width,
          this.options.image.height
        );
      } else {
        // Fall back to the size of the image in the first frame of the first
        // available direction of the default animation if one exists
        const defaultAnimationDirections = Object.values(
          this.options.animations[this.options.defaultAnimation]
        )[0];
        if (
          defaultAnimationDirections &&
          (defaultAnimationDirections.images?.length ?? 0) > 0
        ) {
          this.size = vec(
            defaultAnimationDirections.images![0].width,
            defaultAnimationDirections.images![0].height
          );
        }
      }

      // Otherwise leave the size as (0, 0)
    }

    if (this.options.origin) {
      this.origin = vec.cpy(this.options.origin);
    } else {
      // Default to the center of the sprite based on size
      this.origin = vec.mul(this.size, 0.5);
    }

    if (this.options.scale) {
      this.scale = this.options.scale;
    }

    if (this.options.rotation) {
      this.rotation = this.options.rotation;
    }

    // Check and initialise direction
    this._direction = this.options.defaultDirection;
    if (
      this.options.directions.length === 0 ||
      !this.options.directions.includes(this._direction)
    ) {
      throw new Error(`Invalid direction "${this._direction}"`);
    }

    // Check and initialise animation
    this._animation = this.options.defaultAnimation;
    const animations = Object.keys(this.options.animations);
    if (
      animations.length === 0 ||
      !animations.includes(this._animation)
    ) {
      throw new Error(`Invalid animation "${this._animation}"`);
    }

    // Make sure attachment point keyframes are defined in ascending
    // frame order in all animations
    for (const animation of Object.keys(this.options.animations)) {
      for (const direction of Object.keys(this.options.animations[animation])) {
        if (
          this.options.animations[animation][direction].attachmentPointKeyframes
        ) {
          for (const attachmentPoint of Object.keys(
            this
              .options
              .animations[animation][direction]
              .attachmentPointKeyframes!
          )) {
            this
              .options
              .animations[animation][direction]
              .attachmentPointKeyframes![attachmentPoint]
              .sort(
                (a, b) => a.frame - b.frame
              );
          }
        }
      }
    }
  }

  public get direction(): string {
    return this._direction;
  }

  public set direction(value: string) {
    if (this.options.directions.includes(value)) {
      this._direction = value;
    }
  }

  public get animation(): string {
    return this._animation;
  }

  public set animation(value: string) {
    if (Object.keys(this.options.animations).includes(value)) {
      this._animation = value;
    }
  }

  public playAnimation() {
    if (this.currentAnimationState) {
      this.currentAnimationState.playing = true;
    }
  }

  public pauseAnimation() {
    if (this.currentAnimationState) {
      this.currentAnimationState.playing = false;
    }
  }

  public resetAnimation() {
    if (this.currentAnimationState) {
      this.currentAnimationState.currentFrame = 0;
      this.currentAnimationState.currentFrameTime = 0;
    }
  }

  public getAttachmentPoint(name: string): vec | null {
    return this.currentAttachmentPoints?.[name] ?? null;
  }

  public update(dt: number) {
    this.currentAnimationOptions = this.updateAnimationOptions();
    this.currentAnimationState = this.updateAnimationState(dt);
    this.currentImage = this.updateImage();
    this.currentAttachmentPoints = this.updateAttachmentPoints();
  }

  private updateAnimationOptions(): SpriteAnimationOptions {
    if (!(this._animation in this.options.animations)) {
      throw new Error(`Invalid animation "${this._animation}"`);
    }

    const directions = Object.keys(this.options.animations[this._animation]);
    if (directions.length === 0) {
      throw new Error(
        `No directions available for animation "${this._animation}"`
      );
    }

    if (this._direction in this.options.animations[this._animation]) {
      return this.options.animations[this._animation][this._direction];
    }

    if ('*' in this.options.animations[this._animation]) {
      return this.options.animations[this._animation]['*'];
    }

    return this.options.animations[this._animation][directions[0]];
  }

  private updateAnimationState(dt: number): SpriteAnimationState {
    if (
      !this.currentAnimationOptions ||
      !this.currentAnimationState
    ) {
      return {
        playing: true,
        currentFrame: 0,
        currentFrameTime: 0,
      };
    }

    if (this.currentAnimationState.playing) {
      const frameTime = 1 / this.currentAnimationOptions.frameRate;
      this.currentAnimationState.currentFrameTime += dt;

      if (this.currentAnimationState.currentFrameTime > frameTime) {
        const frameCount = this.currentAnimationOptions.frameCount;
        this.currentAnimationState.currentFrame++;
        this.currentAnimationState.currentFrameTime = 0;

        if (this.currentAnimationState.currentFrame > frameCount) {
          switch (this.currentAnimationOptions.mode) {
            case SpriteAnimationRepeatMode.PlayOnceAndReset:
              this.currentAnimationState.playing = false;
              this.currentAnimationState.currentFrame = 0;
              break;

            case SpriteAnimationRepeatMode.PlayOnceAndStop:
              this.currentAnimationState.playing = false;
              this.currentAnimationState.currentFrame = frameCount - 1;
              break;

            case SpriteAnimationRepeatMode.Repeat:
              this.currentAnimationState.currentFrame = 0;
              break;
          }
        }
      }
    }

    return this.currentAnimationState;
  }

  private updateImage(): HTMLImageElement | HTMLCanvasElement | null {
    if (
      !this.currentAnimationOptions ||
      !this.currentAnimationState
    ) {
      return null;
    }

    if (
      !this.currentAnimationOptions.images ||
      this.currentAnimationOptions.images.length === 0
    ) {
      return this.options.image ?? null;
    }

    return this.currentAnimationOptions.images[
      this.currentAnimationState.currentFrame
    ] ?? this.options.image ?? null;
  }

  private updateAttachmentPoints(): SpriteAttachmentPointMap | null {
    if (
      !this.options.attachmentPoints ||
      this.options.attachmentPoints.length === 0
    ) {
      return null;
    }

    if (!this.currentAttachmentPoints) {
      this.currentAttachmentPoints = Object.fromEntries(
        this.options.attachmentPoints.map(attachmentPoint => [
          attachmentPoint.name,
          attachmentPoint.offset
        ])
      );
    }

    if (
      this.currentAnimationOptions &&
      this.currentAnimationOptions.attachmentPointKeyframes &&
      this.currentAnimationState
    ) {
      for (const name of Object.keys(this.currentAttachmentPoints)) {
        if (
          name in this.currentAnimationOptions.attachmentPointKeyframes &&
          this.currentAnimationOptions.attachmentPointKeyframes[name].length > 0
        ) {
          const previousKeyframe = this.findPreviousKeyframe(
            this.currentAnimationOptions.attachmentPointKeyframes[name],
            this.currentAnimationState.currentFrame
          );
          this.currentAttachmentPoints[name] = previousKeyframe.offset;
        }
      }
    }

    return this.currentAttachmentPoints;
  }

  private findPreviousKeyframe(
    keyframes: SpriteAttachmentPointKeyframe[],
    currentFrame: number
  ): SpriteAttachmentPointKeyframe {
    const found = [...keyframes].reverse().find(
      keyframe => keyframe.frame <= currentFrame
    );

    if (!found) {
      return keyframes[keyframes.length - 1];
    }

    return found;
  }

  public draw(context: CanvasRenderingContext2D) {
    context.save();
    context.translate(
      this.position.x,
      this.position.y
    );
    context.scale(this.scale, this.scale);
    context.rotate(this.rotation);

    this.options.preRender?.(context, this);

    if (this.currentImage) {
      context.drawImage(
        this.currentImage,
        -this.origin.x,
        -this.origin.y,
        this.currentImage.width,
        this.currentImage.height
      );
    }

    this.options.postRender?.(context, this);

    if (this.options.debug.showSpriteBoundingBox) {
      context.strokeStyle = Sprite.DEBUG_BOUNDING_BOX_COLOUR;
      context.lineWidth = Sprite.DEBUG_BOUNDING_BOX_LINE_WIDTH;
      context.strokeRect(
        -this.origin.x,
        -this.origin.y,
        this.size.x,
        this.size.y
      );
    }

    if (this.options.debug.showSpriteTransforms) {
      this.drawTransformsMarker(
        context,
        vec(),
        Sprite.DEBUG_TRANSFORMS_COLOUR_X,
        Sprite.DEBUG_TRANSFORMS_COLOUR_Y,
        Sprite.DEBUG_TRANSFORMS_LINE_WIDTH,
        Sprite.DEBUG_TRANSFORMS_SIZE
      );
    }

    if (
      this.options.debug.showAttachmentPoints &&
      this.currentAttachmentPoints
    ) {
      for (const attachmentPoint of Object.values(this.currentAttachmentPoints)) {
        this.drawCross(
          context,
          attachmentPoint,
          Sprite.DEBUG_ATTACHMENT_POINT_COLOUR,
          Sprite.DEBUG_ATTACHMENT_POINT_LINE_WIDTH,
          Sprite.DEBUG_ATTACHMENT_POINT_SIZE
        );
      }
    }

    context.restore();
  }

  private drawTransformsMarker(
    context: CanvasRenderingContext2D,
    position: vec,
    xColour: string,
    yColour: string,
    lineWidth: number,
    size: number
  ) {
    context.save();

    context.lineWidth = lineWidth;

    context.strokeStyle = xColour;
    context.beginPath();
    context.moveTo(position.x, position.y);
    context.lineTo(position.x + size, position.y);
    context.stroke();

    context.strokeStyle = yColour;
    context.beginPath();
    context.moveTo(position.x, position.y);
    context.lineTo(position.x, position.y + size);
    context.stroke();

    context.restore();
  }

  private drawCross(
    context: CanvasRenderingContext2D,
    position: vec,
    colour: string,
    lineWidth: number,
    size: number
  ) {
    context.save();

    context.lineWidth = lineWidth;

    const halfSize = Math.ceil(size / 2);
    context.strokeStyle = colour;
    context.beginPath();
    context.moveTo(position.x - halfSize, position.y - halfSize);
    context.lineTo(position.x + halfSize, position.y + halfSize);
    context.moveTo(position.x - halfSize, position.y + halfSize);
    context.lineTo(position.x + halfSize, position.y - halfSize);
    context.stroke();

    context.restore();
  }
}

/**
 * Content Manager Processor wrapper which converts SpriteOptionsData into
 * SpriteOptions
 *
 * @see https://www.npmjs.com/package/@basementuniverse/content-manager
 */
export async function spriteOptionsContentProcessor(
  content: Record<string, {
    name: string;
    type: string;
    content: any;
    status: number;
  }>,
  data: {
    name: string;
    type: string;
    content: any;
    status: number;
  }
): Promise<void> {
  const getImageFromContent = (name: string):
    | HTMLImageElement
    | HTMLCanvasElement
    | null => {
    const image = content[name]?.content;
    if (!image) {
      throw new Error(`Image '${name}' not found`);
    }

    return image;
  };

  const result: any = data;
  if (result.imageName) {
    result.image = getImageFromContent(result.imageName);
    delete result.imageName;
  }

  if (result.animations) {
    for (const [animationName, animation] of (
      Object.entries(result.animations) as [
        string,
        {
          [direction: string]: SpriteAnimationOptionsData;
        }
      ][]
    )) {
      for (const [directionName, direction] of Object.entries(animation)) {
        if (direction.imageNames) {
          result.animations[animationName][directionName].images = direction
            .imageNames
            .map(getImageFromContent);
          delete result.animations[animationName][directionName].imageNames;
        }
      }
    }
  }

  data.content = result as SpriteOptions;
}
